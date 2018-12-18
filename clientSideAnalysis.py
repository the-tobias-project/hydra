#!/usr/bin/env python3

# stdlib
import time
import sys
import pdb
import itertools
import _pickle as pickle
import os

# Third party lib
import h5py
import numpy as np
from plinkio import plinkfile
from sklearn.utils.extmath import svd_flip
import tqdm

# Internal lib
from corr import nancorr
from settings import Settings, Commands, Options, QCOptions, PCAOptions, QCFilterNames, PCAFilterNames
from optimizationAux import *
from utils import encode, decode, write_or_replace


class ServerTalker(object):
    """Dispatches commands from/to the server"""
    def __init__(self, store_name, scratch, server):
        self.store_name = store_name
        self.status = ""
        self.scratch = scratch
        self.server = server
        self.store_path = self._store_path()
        if not os.path.exists(scratch):
            os.makedirs(scratch)
        self.local_n = None
        self.chroms = None
        self.do_ld = False
        self.r1 = None
        self.r2 = None
        self.r3, self.r4 = None, None
        self.tqdm = None
        self.verbose = True
        self.covariates = None
        self.Ys = None
        self.previous_estimates = {}
        self.previous_Us = {}

    def _store_path(self):
        plinkName = self.store_name
        basename = os.path.basename(plinkName)
        write_to = self.scratch
        store_name = os.path.join(write_to, basename + '.h5py')
        return store_name

    def plinkToH5(self):
        """Gets plink prefix, produces an HDF file with the same prefix"""
        plinkName = self.store_name
        plink_file = plinkfile.open(plinkName)
        if not plink_file.one_locus_per_row():
            print("""This script requires that snps are
                rows and samples columns.""")
            sys.exit(1)
        sample_list = plink_file.get_samples()
        locus_list = plink_file.get_loci()
        n_tot = len(sample_list)
        with h5py.File(self.store_path, 'w', libver='latest'
                , swmr=True) as store:
            store.attrs['n'] = len(sample_list)
            store.attrs['has_local_AF'] = False
            store.attrs['has_global_AF'] = False
            store.attrs['has_centering'] = False
            store.attrs['has_normalization'] = False
            affection = [sample.affection for sample in sample_list]
            write_or_replace(store, 'meta/Status', affection, np.int8)
            ids = [sample.iid.encode('utf8') for sample in sample_list]
            write_or_replace(store, 'meta/id', ids, 'S11')
            del ids, affection
            # Read Demographic file
            with open(plinkName + ".ind", 'r') as dem_f: 
                dem = [(row.split(",")[2]).encode("UTF8") for row in dem_f]
                write_or_replace(store, 'meta/regions', dem)
            # Read chromosome data
            current_chr = 1
            positions = []
            current_group = store.require_group(str(current_chr))
            genotypes = np.zeros(n_tot, dtype=np.float32)
            for locus in locus_list:
                row = next(plink_file)
                if locus.chromosome != current_chr: 
                    if len(positions) == 0:
                        del store[str(current_chr)]
                    else:
                        write_or_replace(current_group, 'positions', positions,
                            dtype=np.uint)
                        msg = encode({"TASK": "INIT", "SUBTASK": "POS",
                            "CHROM": current_chr, "POS": positions})
                        self.server.message(msg)
                        positions = []
                    current_chr = locus.chromosome
                    if current_chr == 23:
                        break
                    current_group = store.require_group(str(current_chr))
                genotypes[:] = np.array(row, dtype=np.float32)
                pos = str(locus.bp_position)
                dset = current_group.require_dataset(pos, (n_tot,),
                    dtype=np.float32)
                positions.append(pos)
                dset.attrs['rsid'] = locus.name
                dset.attrs['snp'] = locus.allele1
                dset.attrs['alt'] = locus.allele2
                dset.attrs['counts'] = np.array([np.sum(genotypes==1), 
                  np.sum(genotypes==2), np.sum(genotypes==3)], dtype=np.int32)
                genotypes[genotypes==3] = np.nan
                dset[:] = genotypes
            if locus.chromosome != 23:
                write_or_replace(current_group, 'positions', positions,
                    np.uint32)
                msg = {"TASK": "INIT", "SUBTASK": "POS",
                    "CHROM": current_chr, "POS": positions}
                self.server.message(encode(msg))

    def dispatch(self, message):
        task = message["TASK"]
        subtask = message["SUBTASK"]
        if self.verbose:
            print("Working on:", task)
        self.status = "{}/{}".format(task, subtask) 
        if task == Commands.INIT:
            if subtask == "STORE":
                self.plinkToH5()
                print("preparing counts")
                self.reportCounts()
            if subtask == "STATS":
                self.store_stats(message)
        elif task == Commands.QC:
            self.run_QC(message)
        elif task == Commands.PCA:
            # First we will perform all the filters save for LD pruning: 
            if "FILTERS" in message: 
                self.run_snp_filters(message)
            elif self.do_ld: # Just do LD
                self.run_ld(message)
            elif subtask == "PCA_POS":
                if self.chroms is None:
                    self.load_chroms()
                chrom = self.store_message(message, "PCA_passed")
                self.chroms.remove(chrom)
                if not self.chroms:
                    self.chroms = None
                    self.load_chroms()
                    self.standardize_data(scale=False, center=True)
                    self.report_covariance(self.chroms , "PCA_passed")
                    self.r0, self.r1 = None, None
            if subtask == "PCS":
                self.compute_Us(message)
        elif task == Commands.ASSO:
            if self.covariates is None:
                self.load_covariates(message)
            self.run_logistic_regression(message)

    def standardize_data(self, scale=False, center=True):
        with h5py.File(self.store_path, 'a') as store: 
            for chrom in self.chroms:
                group = store[chrom]
                if center:
                    af = group["MAF"].value
                if scale:
                    sd = np.sqrt(group["VAR"].value)
                pos  = group["positions"].value
                for i, snp in enumerate(pos):
                    genotypes = group[str(snp)]
                    if center:
                        vals = genotypes.value
                        vals[np.isnan(vals)] = 2*af[i]
                        genotypes[:] = vals - 2*af[i]
                    if scale:
                        genotypes[:] = genotypes[:] / sd[i]
                        store.attrs["has_normalization"] = True

    def load_chroms(self):
        with h5py.File(self.store_path, 'r') as store:
            self.chroms = [i for i in store.keys() if i!= 'meta']

    def store_stats(self, message):
        chrom = message["CHROM"]
        with h5py.File(self.store_path, 'a') as store: 
            chrom_group = store[str(chrom)]
            if "MISS" in message: 
                vals = message["MISS"]
                task = "not_missing_per_snp"
                dset = chrom_group.create_dataset(task, data=1-vals)
            if "AF" in message:
                vals = message["AF"]
                task = 'MAF'
                dset = chrom_group.create_dataset(task, data=vals)
            if "HWE" in message:
                vals = message["HWE"]
                task = "hwe"
                dset = chrom_group.create_dataset(task, data=vals)
            if "VAR" in message:
                vals = message["VAR"]
                task = "VAR"
                dset = chrom_group.create_dataset(task, data=vals)

    def reportCounts(self):
        """Report the counts (Het, homo Alt, missing)"""
        with h5py.File(self.store_path, 'r') as store:
            countDict = {}
            n = store.attrs["n"]
            countDict["START"] = True
            keys = [i for i in store.keys() if i != 'meta']
            for chrom in keys:
                countDict["n"] = int(n)
                countDict["TASK"] ="INIT"
                countDict["SUBTASK"] = "COUNT"
                countDict["CHROM"] = chrom
                if chrom == 'meta':
                    continue
                positions = store["{}/positions".format(chrom)].value
                count_arr = []
                for i, pos in enumerate(positions):
                    count_arr.append(store["{}/{}".format(chrom, 
                        pos)].attrs["counts"])
                countDict["COUNTS"] = np.array(count_arr, dtype=np.uint32)
                if chrom == keys[-1]:
                    countDict["END"] = True
                self.server.message(encode(countDict))
                countDict = {}

#################################### QC #######################################
    def run_QC(self, message, remove=True):
        def find_what_passes(qc_name, dset_name, tokeep, doubleSided=False):
            vals = group[dset_name].value
            if qc_name in filter_list:
                ind = filter_list.index(qc_name)
                if not doubleSided:
                    tokeep = np.logical_and(tokeep, vals > value_list[ind])
                else:
                    tokeep = np.logical_and(tokeep,
                        np.logical_and(vals > value_list[ind], 
                            (1.0-vals) > value_list[ind]))
            return tokeep

        def replace_dataset(tokeep, dset_name, return_deleted=False):
            vals = group[dset_name].value
            remaining = vals[tokeep]
            deleted = vals[np.logical_not(tokeep)]
            write_or_replace(group, dset_name, remaining)
            if return_deleted:
                return deleted

        filter_list = message["FILTERS"]
        value_list = message["VALS"]
        with h5py.File(self.store_path, 'a') as store:
            if self.chroms is None:
                self.chroms = [i for i in store.keys() if i != "meta"]

            for chrom in self.chroms:
                group = store[chrom]
                positions = group['positions'].value
                tokeep = np.ones_like(positions, dtype=bool)
                tokeep = find_what_passes(QCFilterNames.QC_HWE, "hwe", tokeep)
                tokeep = find_what_passes(QCFilterNames.QC_MAF, "MAF",
                    tokeep, doubleSided=True)
                if QCFilterNames.QC_MPS in filter_list:
                    ind = filter_list.index(QCFilterNames.QC_MPS)
                    value_list[ind] = 1 - value_list[ind]
                tokeep  = find_what_passes(QCFilterNames.QC_MPS, "not_missing_per_snp", tokeep)
                print("After filtering, {} snps remain".format(np.sum(tokeep)))
                if remove: # Delete what doesn't pass 
                    replace_dataset(tokeep, 'hwe')
                    replace_dataset(tokeep, 'MAF')
                    replace_dataset(tokeep, 'not_missing_per_snp')
                    deleted = replace_dataset(tokeep, 'positions', 
                        return_deleted=True)
                    for snp in deleted:
                        snp = str(snp)
                        if snp in group:
                            del group[snp]
                else: # Store what has been tagged
                    if "passed" in group: 
                        tags = group["passed"]
                        tokeep = np.logical_and(tokeep, tags.value)
                        tags[:] = tokeep
                    else:
                        if "PCA_passed" in group:
                            del group["PCA_passed"]
                        if "PCA_positions" in group:
                            del group["PCA_positions"]
                        group.create_dataset("PCA_passed", 
                            data=np.ones(np.sum(tokeep), dtype=bool))
                        positions = group['positions'].value[tokeep]
                        group.create_dataset("PCA_positions", data=positions)

#################################### PCA ######################################
    def run_snp_filters(self, message):
        print("Preparing for PCA")
        filters = message["FILTERS"]
        vals    = message["VALS"]
        if PCAFilterNames.PCA_LD in filters: 
            ind = filters.index(PCAFilterNames.PCA_LD)
            winsize = vals.pop(ind)
            filters.pop(ind)
            self.do_ld = True # deal with LD after other filters
            self.r1 = int(winsize)
            self.r2 = int(winsize/2)
        self.run_QC({"FILTERS": filters, "VALS": vals}, remove=False)
        if self.do_ld:
            self.r3 = 0
            with h5py.File(self.store_path, 'r') as fp: 
                if self.chroms is None:
                    self.chroms = [i for i in store.keys() if i!= 'meta']
                n = 0
                for chrom in self.chroms:
                    n = max(n, len(fp["{}/PCA_positions".format(chrom)]))

            print("LD pruning. This will take some time...")
            self.tqdm = tqdm.tqdm(total=n)
            self.verbose = False
            self.run_ld({})
 
    def run_ld(self, message):
        msg = {"TASK": "PCA", "SUBTASK": "LD"}
        with h5py.File(self.store_path, 'r') as store:
            n = store.attrs['n']
            for chrom in self.chroms:
                tags = store["{}/PCA_passed".format(chrom)]
                if chrom in message:
                    state = message[chrom]
                    if state[0] == "E": # Finished with this chrom
                        self.chroms.remove(chrom)
                        if len(self.chroms) == 0:
                            self.do_ld = False 
                            self.chroms = None
                            msg = {"TASK": "PCA", "SUBTASK": "PCA_POS"}
                            self.server.message(encode(msg))
                            self.verbose = True
                            return
                        continue
                    else:
                        tokeep = state
                        end = self.r3 + len(tokeep)
                else:
                    end = min(self.r3 + self.r1, tags.shape[0])
                    tokeep = tags[self.r3: end]
                pos = store["{}/PCA_positions".format(chrom)]
                positions = pos[self.r3: end]
                positions = positions[tokeep] 
                genotypes = np.empty((n,len(positions)), dtype=np.float32)
                for i, snp in enumerate(positions): 
                    genotypes[:,i] = store["{}/{}".format(chrom, snp)].value
                corr = nancorr(genotypes)
                msg[chrom] = corr
        self.server.message(encode(msg))
        self.r3 += self.r2
        self.tqdm.update(self.r2)

    def store_message(self, message, dset_name):
        for key, val in message.items():
            if key in ["TASK", "SUBTASK"]:
              continue 
            else: 
              with h5py.File(self.store_path, 'a') as store:
                  dset = store.require_dataset("{}/{}".format(key, dset_name),
                      val.shape, dtype = val.dtype)
                  dset[:] = val
              return key

    def report_covariance(self, chroms, mask_name):
        print("reporting cov")
        msg = {"TASK": "PCA", "SUBTASK": "COV"}
        with h5py.File(self.store_path, 'r') as store:
            n = store.attrs["n"]
            for i_ch1, ch1 in enumerate(chroms):
                group = store[ch1]
                pos   = group["PCA_positions"].value
                mask  = group[mask_name].value
                pos   = pos[mask]
                g1     = np.empty((len(pos), n))

                for i, snp1 in enumerate(pos):
                    g1[i, :] = group[str(snp1)].value

                for i_ch2, ch2 in enumerate(chroms):
                    if i_ch2 > i_ch1:
                        break
                    group = store[ch2]
                    pos   = group["PCA_positions"].value
                    mask  = group[mask_name].value
                    pos   = pos[mask]
                    g2     = np.empty((n, len(pos)))
                    for i, snp2 in enumerate(pos):
                        g2[:, i] = group[str(snp2)].value
                    msg["CH1"] = ch1
                    msg["CH2"] = ch2
                    msg["MAT"] = g1.dot(g2).astype(np.float32)
                    if ch1 == chroms[-1] and ch2 == chroms[-1]:
                        msg["E"] = True
                    self.server.message(encode(msg))

    def compute_Us(self, message):
        inv_sigma = message["ISIG"]
        v = message["V"]
        chroms = message["CHROMS"]
        with h5py.File(self.store_path, 'a') as store:
            dset = store["meta"]
            pca_sigma = dset.require_dataset('pca_sigma', shape=inv_sigma.shape, dtype=np.float32)
            n = 0
            for chrom in chroms: 
                n += np.sum(store["{}/PCA_passed".format(chrom)])
            num_inds = store.attrs["n"]
            arr = np.empty((num_inds, n))
            offset = 0
            for chrom in chroms:
                group = store[str(chrom)]
                tokeep    = group["PCA_passed"].value
                positions = group["PCA_positions"].value[tokeep]
                for i, position in enumerate(positions): 
                    arr[:, offset+i] = group[str(position)].value
                offset += i
            u = arr.dot(v.T).dot(np.diag(inv_sigma))
            u, v = svd_flip(u, v, u_based_decision=False)
            pca_vt = dset.require_dataset('pca_v.T', shape=v.shape,
                dtype=np.float32)
            pca_vt[:,:] = v
            pca_u = dset.require_dataset('pca_u', shape=u.shape,
                dtype=np.float32)
            pca_u[:,:] = u

################################ ASSOCIATION ##################################
    def load_covariates(self, message):
        vals = message["VARS"]
        with h5py.File(self.store_path, 'r') as store:
            n = store.attrs["n"]
            if self.Ys is None:
                self.Ys = np.sign(store["meta/Status"].value 
                    - 0.5).reshape(n, 1)
            self.covariates = np.empty((n,vals[0] + len(vals)))
            if vals[0] > 0: # Load up scores 
                scores = store["meta/pca_u"].value[:,:vals[0]]
                scores /= np.std(scores, axis = 0)
                self.covariates[:,1:vals[0]+1] = scores  # already centered 
                self.covariates[:,1:] *= -self.Ys

    def run_logistic_regression(self, message):
        if "VARS" in message:
            self.ncov = message["VARS"][0] + 1 #Number of Pcs plus the snp 
        if "CHROM" in message:
            chrom = message["CHROM"]
            warm_start = message["VALS"]
            msg = self._run_logistic_regression(chrom, self.ncov, warm_start)
        else: # First iteration
            self.load_chroms()
            for chrom in self.chroms:
                msg = self._run_logistic_regression(chrom, self.ncov,  None)
        self.server.message(msg)

    def _run_logistic_regression(self, chrom, ncov, 
            warm_start=None, rho=10.0, alpha=1.0):
        # ASSUMES EVERYTHING IS CENTERED #TODO 
        covariates = self.covariates
        with h5py.File(self.store_path, 'r') as store:
            group = store[chrom]
            positions = group["positions"]
            estimates = np.empty((ncov, len(positions))) # Can probably get away with float32 here
            estimates[:, -1] = 0
            is_standardized = store.attrs["has_normalization"]
            std = np.sqrt(group["VAR"])
            offset = 1
            if chrom in self.previous_estimates:
                z_hat = self.previous_estimates[chrom]
                all_Us     = self.previous_Us[chrom] + z_hat - warm_start
            else:
                all_Us     = 0
            for i, position in enumerate(positions):
                if std[i] == 0:
                    estimates[:, i] = np.nan
                    offset += 1
                    continue
                if is_standardized:
                    gval = group[str(position)].value
                else:
                    gval = group[str(position)].value/std[i]
                covariates[:,0] = gval
                if warm_start is None:
                    estimates[:,i] = bfgs_more_gutted(covariates,
                        np.zeros((ncov, )), np.zeros((ncov,)),
                        rho, estimates[:,i-offset], ncov)
                    z_hat = alpha * estimates
                else:
                    estimates[:,i] = bfgs_more_gutted(covariates,
                        all_Us[:,i] , warm_start[:,i], rho, z_hat[:, i], ncov)
                    z_hat = alpha * estimates + (1-alpha) * warm_start
                offset = 1
        self.previous_estimates[chrom] = estimates
        self.previous_Us[chrom] = all_Us
        return encode({"TASK": Commands.ASSO, "SUBTASK": chrom, 
          "CHROM": chrom, "VALS": z_hat})


if __name__=='__main__':
    print("no commands here yet. Test using WTCCC_run.py")
