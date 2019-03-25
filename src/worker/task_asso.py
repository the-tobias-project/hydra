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

class LogisticAdmm(object):
    __instance = None
    def __init__(self, covar_numbers, npcs):
        if LdReporter.__instance is not None:
            return 
        else:
            LdReporter.__instance = self
            pfile = client_config["plinkfile"]
            store_name = shared.get_plink_store(pfile)
            self.store = h5py.File(store_name, 'a')
            self.load_Y(covar_numbers[-1])
            self.load_covar_mat(npcs, covar_numbers[:-1])
            LdReporter.__instance = self

    @staticmethod
    def get_instance(covar_numbers, npcs):
        if LdReporter.__instance is None:
            LdReporter(covar_numbers, npcs)
        return LdReporter.__instance
    
    def load_Y(self, location):
        name = shared.get_covar_file()
        y = np.loadtxt(name, usecols=[location], delimiter='\t')
        n = y.shape[0]
        self.Ys = np.sign(y - .5).reshape(n, 1)

    def load_covar_mat(self, npcs, other_covars_loc):
        store = self.store
        n = store.attrs["n"]
        p = npcs + len(other_covars) + 2 # intercept + snp
        self.covariates = np.empty((n, p))
        # Load up scores 
        scores = store["meta/pca_u"].value[:,:npcs]
        scores /= np.std(scores, axis = 0)  # This needs to be the global std
        self.covariates[:, 2:npcs+1] = -self.Ys * scores  # already centered
        del scores
        covar_file_name = shared.get_covar_file()
        other_covars = np.loadtxt(covar_file_name, delimiter="\t", usecols=other_covars_loc)
        other_covars /= np.std(other_covars, axis = 0)  # This needs to be the global std
        other_covars -= np.mean(other_covars, axis = 0)  # This needs to be the global std
        self.covariates[:, npcs+1:] = -self.Ys * other_covars
        self.covariates[:, 0] = -self.Ys.reshape(-1)

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
