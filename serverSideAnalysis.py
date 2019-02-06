#!/usr/bin/env python 3

# stdlib
import sys
import os
import time
import _pickle as pickle
import logging
import pdb

# Third party lib
import tqdm
import numpy as np
import h5py
from plinkio import plinkfile
from scipy.sparse.linalg import eigsh as eig

# Internal lib
from corr import corr, hweP
from settings import Settings, Commands, Options, QCOptions, PCAOptions, QCFilterNames, PCAFilterNames
from utils import encode, decode, write_or_replace

#TODO documentation
# Careful here, eigh uses https://software.intel.com/en-us/mkl-developer-reference-c-syevr behind the hood
# so it can be significantly slower 


class ServerTalker(object):
    """Dispatches and analyses of commands from the server."""
    def __init__(self, connections, server, scratch=Settings.local_scratch):
        self.connections = connections
        self.status      = "Initialized"
        self.task        = None
        self.subtasks    = None
        self.server      = server
        self.scratch     = scratch
        self.storePath   = os.path.join(scratch, "central.h5py")
        self.store       = h5py.File(self.storePath, "a")
        self.counter     = connections
        self.r0          = None # General, multipurpose registers
        self.r1, self.r2 = None, None
        self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
        self.tqdm        = None
        self.estimates   = {}
        self.max_iters   = 5
        self.finished    = set()
        self.iters       = {}
        self.logger      = logging.getLogger(__name__)

    def report_status(self):
        self.logger.info("Now working on: {}".format(self.status))

    def dispatch(self,message):
        task = message["TASK"]
        self.task = task
        subtask = message["SUBTASK"]
        self.subtask = subtask
        current_task = "{} {}".format(task, subtask)
        self.logger.debug(current_task)
        if self.status != current_task:
            self.status = current_task
            self.report_status()
        if self.task == Commands.INIT:
            if subtask == "POS":
                self.store_positions(message)
            elif subtask == "COUNT":
                self.store_counts(message)
                self.logger.info("Count statistics have been initialized!")
        elif self.task == Commands.QC:
            self.QC_filters(message["SUBTASK"], message["VALS"])
        elif self.task == Commands.PCA:
            if self.subtask == "FILTERS":
                continue_to_ld_prune = self.PCA_filters(message["FILTERS"]
                    , message["VALS"])
                self.counter = self.connections
                self.r0      = 0
                if continue_to_ld_prune:
                    return
            if self.subtask == "LD":
                self.ld_filters(message)
            if self.subtask == "PCA_POS":
                self.counter -= 1
                if self.counter == 0: #NOTE I don't think this is neccessary
                    msg = {"TASK": self.task, "SUBTASK": self.subtask}
                    self.update_pca_positions()
                    self.report_content("PCA_passed", msg)
                    self.counter = None
            if self.subtask == "COV":
                self.logger.info("building covariance")
                if self.buildCov(message):
                  return
      #              self.pca()
      #              self.server.get_options()
      #  elif self.task == Commands.ASSO:
      #      self.run_logistic_regression(message)

    def store_positions(self, message):
        chrom = message["CHROM"]
        positions = message["POS"]
        dsetname = "{}/positions".format(chrom)
        write_or_replace(self.store, dsetname, positions, np.uint32)
        self.logger.info("{} loci in chromosome {}.".format(len(positions), chrom))

    def store_counts(self, message):
        n = message["n"]
        if "START" in message:
            if "N" not in self.store.attrs:
                self.store.attrs["N"] = 0
            self.store.attrs["N"] += n
        chrom = message["CHROM"]
        size = len(message["COUNTS"])
        dsetname = "{}/counts".format(chrom)
        if dsetname not in self.store:
            dset = self.store.require_dataset(dsetname, (size, 4)
              , dtype=np.uint32)
        else:
            dset = self.store[dsetname]
        counts = message["COUNTS"]
        homo_ref = n - np.sum(counts, axis = 1)[:, np.newaxis].astype(
            np.uint32)
        dset[:] += np.hstack((homo_ref, counts))
        if "END" in message:
            self.counter -= 1
            if self.counter == 0: # Everybody's report is in
                self.status = "INIT STATS"
                self.report_status()
                self.count_stats()
                self.server.get_options()


    def count_stats(self):
        N = float(self.store.attrs["N"])
        task = "INIT"
        for chrom in self.store.keys():
            counts_dset = self.store["{}/counts".format(chrom)].value
            missing_rate = counts_dset[:,3] / float(N)
            missing_rate_dset = self.store.create_dataset(
                "{}/missing_rates".format(chrom), data=missing_rate)
            af = (counts_dset[:,2] * 2 + counts_dset[:,1]).astype(float)
            af /= (np.sum(counts_dset[:,:3], axis=1)*2).astype(float)
            #af = np.minimum(af, 1-af)
            self.store.create_dataset("{}/allele_freq".format(chrom), data=af)
            #var = counts_dset[:,0] * (2*af)**2 
            #var += counts_dset[:,1] * (1-2*af)**2 
            #var += counts_dset[:,2] * (2-2*af)**2
            #var /= (N-counts_dset[:,3]) # 2*af*(1-af)
            var = 2*af*(1-af)
            self.store.create_dataset("{}/var".format(chrom), data=var)
            hwe = hweP(counts_dset[:,:3].astype(np.int32), 1, 0)
            # Need to Recompile HWEP with uint32
            time.sleep(1)
            msg = {"TASK": task, "SUBTASK": "STATS", "CHROM": chrom
                , "HWE": hwe, "MISS": missing_rate, "AF": af, "VAR": var}
            #pdb.set_trace()
            self.server.message(encode(msg))
            hwe_dset = self.store.create_dataset("{}/hwe".format(chrom)
                , data=hwe)
        time.sleep(1)

################################### QC #######################################
    def QC_filters(self, filter_list, value_list, prefix=None):
        for chrom in self.store.keys():
            if chrom == 'meta':
                continue
            group  = self.store[chrom]
            pos    = group['positions']
            counts = group['counts']
            mr     = group['missing_rates']
            af     = group['allele_freq']
            hwe    = group['hwe']
            tokeep = np.ones(shape=pos.value.shape, dtype=bool)
            if QCFilterNames.QC_HWE in filter_list:
                ind = filter_list.index(QCFilterNames.QC_HWE)
                tokeep = np.logical_and(tokeep, hwe.value > value_list[ind])
            if QCFilterNames.QC_MAF in filter_list:
                ind = filter_list.index(QCFilterNames.QC_MAF)
                tokeep = np.logical_and(tokeep, af.value > value_list[ind] - Settings.kSmallEpsilon)
                tokeep = np.logical_and(tokeep, 1.0-af.value > value_list[ind] - Settings.kSmallEpsilon)
            if QCFilterNames.QC_MPS in filter_list:
                ind = filter_list.index(QCFilterNames.QC_MPS)
                tokeep = np.logical_and(tokeep, mr.value < value_list[ind])
            self.logger.info("in chromosome {}, {} snps were deleted and {} snps remain".format(chrom,
                tokeep.shape[0] - np.sum(tokeep), np.sum(tokeep)))
            # Delete or tag the filtered locations
            if prefix is None: 
                pos_vals, counts_vals = pos.value[tokeep], counts.value[tokeep]
                mr_vals = mr.value[tokeep]
                del group["positions"], group["counts"], group["missing_rates"]
                d1 = group.require_dataset("positions", pos_vals.shape
                    , dtype = pos_vals.dtype)
                d1[:] = pos_vals
                d2 = group.require_dataset("counts", counts_vals.shape
                    , dtype=counts_vals.dtype)
                d2[:] = counts_vals
                d3 = group.require_dataset("missing_rates", mr_vals.shape
                    , dtype=mr_vals.dtype)
                d3[:] = mr_vals
                del pos_vals, counts_vals, mr_vals
                af_vals, hwe_vals= af.value[tokeep], hwe.value[tokeep]
                del group["hwe"], group["allele_freq"]
                d4 = group.require_dataset("hwe", hwe_vals.shape
                    , dtype=hwe_vals.dtype)
                d4[:] = hwe_vals
                d5 = group.require_dataset("allele_freq", af_vals.shape
                    , dtype=af_vals.dtype)
                d5[:] = af_vals
            else:
                n = np.sum(tokeep)
                ones = np.ones(n, dtype=bool)
                d1 = group.require_dataset(prefix + "passed", ones.shape
                    , dtype=bool)
                d1[:] = ones
                pos_vals = pos.value[tokeep]
                d2 = group.require_dataset(prefix + "positions", pos_vals.shape
                    , dtype=pos_vals.dtype)
                d2[:] = pos_vals
                af_vals = af.value[tokeep]
                d3 = group.require_dataset(prefix + "allele_freq"
                    , af_vals.shape, dtype=af_vals.dtype)
                d3[:] = af_vals

################################### PCA #######################################

    def PCA_filters(self, filter_list, value_list):
        ld_prune = False
        if PCAFilterNames.PCA_NONE in filter_list:
            return ld_prune
        if PCAFilterNames.PCA_LD in filter_list: 
            ld_prune = True
            ind = filter_list.index(PCAFilterNames.PCA_LD)
            winsize = value_list.pop(ind)
            filter_list.pop(ind)
            self.r1 = int(winsize)
            self.r2 = int(winsize/2)
        if not set(filter_list).issubset(set(QCOptions.all_options)):
            raise ValueError("At least one of the filters provided is not supported")
        self.QC_filters(filter_list, value_list, prefix="PCA_")
        if ld_prune:
            n = 0
            chroms = [v for v in self.store.keys() if v != 'meta']
            for chrom in chroms:
                n = max(n, len(self.store["{}/PCA_positions".format(chrom)]))
            self.tqdm = tqdm.tqdm(total=n)
        return ld_prune

    def ld_filters(self, message, thresh=0.2):
        self.counter -= 1
        msg = {"TASK":"PCA", "SUBTASK":"LD"}
        for key, val in message.items():
            if key == "TASK" or key == "SUBTASK":
                continue
            else:
                if key in self.sumLin:
                    self.sumLin[key] += val[0]
                    self.sumSq[key]  += val[1]
                    self.cross[key]  += val[2]
                else:
                    self.sumLin[key] = val[0]
                    self.sumSq[key]  = val[1]
                    self.cross[key]  = val[2]
        if self.counter == 0: # Every report is in
            for key, val in message.items():
                if key == "TASK" or key == "SUBTASK":
                    continue
                corr_tot = corr(self.sumLin[key], self.sumSq[key],
                    self.cross[key])
                group  = self.store[key]
                tokeep = self.store["{}/PCA_passed".format(key)].value
                #positions = self.store["{}/PCA_positions".format(key)].value
                end = min(self.r1 + self.r0, len(tokeep))
                maf = self.store["{}/PCA_allele_freq".format(key)].value[
                    self.r0:end]
                maf = maf[tokeep[self.r0:end]]
                maf = 1-np.minimum(maf, 1-maf)
                n = maf.shape[0]
                unfiltered = np.ones((n, ), dtype=bool)
                while True:
                    #length_of_window = np.sum(tokeep[self.r0:end)
                    for i, snp1 in enumerate(unfiltered):
                        if not snp1: # already filtered
                            continue
                        else:
                            for j in range(i+1,n):
                                snp2 = unfiltered[j]
                                if not snp2: # if it didn't pass the filters
                                    continue
                                elif corr_tot[i,j] ** 2 > thresh:
                                    if maf[i] > (maf[j] * (1.0 + 
                                        Settings.kSmallEpsilon)):
                                        unfiltered[i] = False
                                    else:
                                        unfiltered[j] = False
                                    break
                    r2 = corr_tot[unfiltered,:][:,unfiltered]
                    if np.max(r2**2) < thresh:
                        break
                tokeep[self.r0:end][tokeep[self.r0:end]] = unfiltered
                pca_passed = self.store["{}/PCA_passed".format(key)]
                pca_passed[:] = tokeep
                end = min(self.r1 + self.r0 + self.r2, len(tokeep))
                if self.r2 + self.r0 >= len(tokeep) - 1:
                    msg[key] = "END"
                else:
                    msg[key] = tokeep[self.r0 + self.r2:end]
            self.server.message(encode(msg))
            self.counter = self.connections
            self.r0 += self.r2
            self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
            self.tqdm.update(self.r2)

    def update_pca_positions(self):
        chroms = [v for v in self.store.keys() if v != 'meta']
        for chrom in chroms:
            group = self.store[chrom]
            pca_pos = group["PCA_positions"].value
            pca_passed = group["PCA_passed"].value
            write_or_replace(group, "PCA_positions", pca_pos[pca_passed])

    def report_content(self, dset_name, msg):
        chroms = [v for v in self.store.keys() if v != 'meta']
        original_msg = msg.copy()
        for chrom in chroms: 
            dset = self.store["{}/{}".format(chrom, dset_name)]
            msg[chrom] = dset.value
            self.server.message(encode(msg))
            time.sleep(1)
            msg = original_msg.copy()

    def buildCov(self, msg):
        # store the chunks in the store. build on them 
        if self.counter is None:
            self.counter = self.connections
        ch1 = msg["CH1"] 
        ch2 = msg["CH2"]
        if "meta" not in self.store:
            self.store.create_group("meta")
        group = self.store["meta"]
        data = msg["MAT"]
        cov_name = "{}_{}".format(ch1, ch2)
        if cov_name in group:
            data += group[cov_name].value
        write_or_replace(group, cov_name, data)
        print(cov_name)
        if "E" in msg:
            self.counter -= 1
        if self.counter == 0:
            self.logger.info("All covariances have been reported")
            return True
        return False

    def pca(self, chroms=None, n_components=4):
        if chroms is None: 
            chroms = sorted([v for v in self.store.keys() if v != 'meta'])
        cov_size = 0
        meta = self.store["meta"]
        for chrom in chroms: 
            cov_size += meta["{}_{}".format(chrom, chrom)].shape[0]
        self.logger.info("Starting covariance matrix of size {} x {}".format(
            cov_size, cov_size))
        cov = np.empty((cov_size, cov_size))
        i_old = 0
        for chrom1 in chroms:
            j_old = 0
            for chrom2 in chroms:
                if chrom2 > chrom1: 
                    break
                cov_name = "{}_{}".format(chrom1, chrom2)
                if cov_name in meta:
                    pcov = meta[cov_name].value
                cov[i_old:i_old+pcov.shape[0], j_old:j_old+pcov.shape[1]] = pcov
                cov[j_old:j_old+pcov.shape[1], i_old:i_old+pcov.shape[0]] = pcov.T
                j_old += pcov.shape[1]
            i_old += pcov.shape[0]
        cov /= (cov.shape[0])
        sigma, v = eig(cov, k=n_components, ncv=3*n_components)
        sigma, v = zip(*sorted(zip(sigma, v.T), reverse=True))
        v = np.array(v)
        sigma = np.array(sigma)
        sigma[sigma < 0] = 0
        self.store['meta'].create_dataset('Sigmas', data = sigma)
        self.store['meta'].create_dataset('Vs', data = v)
        sigma = np.sqrt(sigma) * np.sqrt(cov.shape[0])
        inv_sigma = sigma.copy()
        inv_sigma[inv_sigma>0] = 1 / inv_sigma[inv_sigma > 0]
        msg = {"TASK": "PCA", "SUBTASK": "PCS", "ISIG": inv_sigma, "V": v
            , "CHROMS": chroms}
        self.server.message(encode(msg))
        time.sleep(0.1)

################################ ASSOCIATION ###################################
    def run_logistic_regression(self, message):
        chrom = message["CHROM"]
        if chrom in self.finished:
            return
        z_hat = message["VALS"]
        if chrom in self.estimates:
            prev = self.estimates[chrom]
            if prev[1] == 1: 
                beta = (prev[0] + z_hat)/(self.connections)
                self.iters[chrom] += 1
                if self.iters[chrom] == self.max_iters:
                    write_or_replace(self.store, chrom + "/results", beta)
                    del self.estimates[chrom]
                    self.finished.add(chrom)
                    chroms = [key for key in self.store if key != 'meta']
                    if len(self.finished) == len(chroms):
                        self.logger.info("We are all done with the regression")
                        self.server.get_response(Commands.all_commands)
                self.estimates[chrom] = [beta, self.connections - 1]
                self.logger.info("Finished iteration{1} on chrom {}".format(self.iters[chrom], chrom))
                self.server.message(encode({"TASK": Commands.ASSO, 
                  "SUBTASK": None, "CHROM": chrom, "VALS": beta}))
            else:
                self.estimates[chrom] = [prev[0] + z_hat , prev[1]-1]

        else: # Not in dictionary yet
            self.estimates[chrom] = [z_hat, self.connections - 1]
            self.iters[chrom] = 1

if __name__=='__main__':
   print("no commands here yet. Test using WTCCC_run.py")
