# stdlib
import os
import pickle
import socket
import sys
import time

# third party lib
from celery import current_app
import h5py
import numpy as np
from plinkio import plinkfile


# internal lib
from client.lib import shared
from lib import HTTPResponse
from lib.utils import write_or_replace
from lib.corr import nancorr, process_plink_row

class LdReporter(object):
    __instance = None

    def __init__(self, win_size, client_config):
        if LdReporter.__instance is not None:
            return 
        else:
            self.r3 = 0
            self.r1, self.r2 = win_size, int(win_size/2)
            LdReporter.__instance = self
            pfile = client_config["plinkfile"]
            self.store = h5py.File(shared.get_plink_store(pfile), 'a')
            self.chroms = [key for key, items in self.store.items() if key != "meta"]
            LdReporter.__instance = self

    @staticmethod
    def get_instance(win_size, client_config):
        if LdReporter.__instance is None:
            return LdReporter(win_size, client_config)
        return LdReporter.__instance

    def update(self, data, client_config):
        data = pickle.loads(data)
        n = self.store.attrs['n']
        msg = {}
        if self.r3 == 0: # fake start the first round
            for chrom in self.chroms: 
                # if the length is less than r1, you deserve an error. 
                # No apologies 
                tags = self.store["{}/PCA_passed".format(chrom)]
                data[chrom] = tags[0:self.r1]
        for key, state in data.items():
            if key == "TASK" or key == "SUBTASK":
                continue
            chrom = key
            tags = self.store["{}/PCA_passed".format(chrom)]
            if state[0] == "E": # Finished with this chrom
                if len(data) == 1: # Done with everything
                    msg = pickle.dumps({})
                    HTTPResponse.respond_to_server('api/tasks/PCA/PCAPOS', 'POST', msg,
                        client_config['name'])
                    print("Done with LD pruning")
                    return
                continue
            else:
                tokeep = state
                end = self.r3 + len(tokeep)
            pos = self.store["{}/PCA_positions".format(chrom)]
            positions = pos[self.r3: end]
            positions = positions[tokeep]
            genotypes = np.empty((n,len(positions)), dtype=np.float32)
            for i, snp in enumerate(positions):
                genotypes[:,i] = self.store["{}/{}".format(chrom, snp)].value
            corr = nancorr(genotypes)
            msg[chrom] = corr
        msg = pickle.dumps(msg)
        HTTPResponse.respond_to_server('api/tasks/PCA/LD', 'POST', msg, client_config['name'])
        self.r3 += self.r2
        print(self.r3)


def store_filtered(message, client_config):
    pfile = client_config["plinkfile"]
    msg = pickle.loads(message)
    with h5py.File(shared.get_plink_store(pfile), 'a') as store:
        for chrom, val in msg:
            dset = store.require_dataset(f"{key}/PCA_passed", val.shape, dtype = val.dtype)
            dset[:] = val


def report_cov(client_config):
    def standardize_mat(mat, af, sd):
        af = 2 * af.reshape(af.shape[0], 1)
        mat -= af
        ind = sd>0
        mat[ind, :] /= sd[ind].reshape(np.sum(ind),1)
        mat[np.isnan(mat)] = 0
        return mat
    pfile = client_config["plinkfile"]
    msg = {}
    with h5py.File(pfile, 'r') as store:
        n = store.attrs["n"]
        for ch1 in store:
            if ch1 == "meta":
                continue
            group = store[ch1]
            tokeep = group['PCA_mask'].value
            pos = group["positions"].value[tokeep]
            af1 = group["MAF"].value[tokeep]
            sd1 = np.sqrt(group["VAR"].value[tokeep])
            g1 = np.empty((len(pos), n))
            for i, snp1 in enumerate(pos):
                g1[i, :] = group[str(snp1)].value
            g1 = standardize_mat(g1, af1, sd1)
            for ch2 in store: 
                if ch2 == "meta":
                    continue
                print(f"Reporting cov: {ch1}_{ch2}")
                group = store[ch2]
                tokeep = group['PCA_mask'].value
                af2 = group["MAF"].value[tokeep]
                sd2 = np.sqrt(group["VAR"].value[tokeep])
                pos   = group["positions"].value[tokeep]
                g2     = np.empty((n, len(pos)))
                for i, snp2 in enumerate(pos):
                    g2[:, i] = group[str(snp2)].value
                g2 = standardize_mat(g2.transpose(), af2, sd2).transpose()
                msg["CH1"] = ch1
                msg["CH2"] = ch2
                msg["MAT"] = g1.dot(g2).astype(np.float32)
                if "E" in message:
                    msg["E"] = True
                print(ch1, ch2)
                msg = pickle.dumps(msg)
                HTTPResponse.respond_to_server('api/tasks/PCA/COV', 'POST', msg, client_config['name'])
