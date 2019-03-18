# stdlib
import logging
import pickle
import re
import os
import pdb

# third party lib
import requests
import h5py
import numpy as np

# internal lib
from lib.settings import Settings, Options, PCAFilterNames, Commands
from lib.client_registry import Registry
from . import task_qc
from lib.corr import corr


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


def decomposed():
    return ("meta" in store and "Sigmas" in store["meta"])

def filtered():
    chrom = [key for key in store if key != "meta"][0]
    if "PCA_passed" in store[chrom] :
        return True
    return False

def start_pca_filters(filters):
   logging.info("Initiating PCA filters...")
   qc_subset_filters = {key: value for key,value in filters.items() 
      if key != PCAFilterNames.PCA_LD}
   if qc_subset_filters: # Perform the qc like filters 
      task_qc.start_local_qc_task(qc_subset_filters, prefix="PCA_")
      filters["remove"] = False
      task_qc.start_client_qc_task(filters, stage=Commands.PCA)


class CovarianceAggregator(object):
    __instance = None
    def __init__(self, num_clients, win_size=50, thresh=0.2):
        if CovarianceAggregator.__instance is not None:
            return 
        else:
            self.num_clients = num_clients
            self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
            self.counter = 0
            self.r0 = 0
            self.r1 = int(win_size)
            self.r2 = int(win_size/2)
            self.send_request({})
            CovarianceAggregator.__instance = self
            self.thresh = thresh

    @staticmethod
    def get_instance(num_clients, win_size):##TODO this is dangerous. If num_clients or win_size changes
        if CovarianceAggregator.__instance is None:
            return CovarianceAggregator(num_clients, win_size)
        return CovarianceAggregator.__instance
    
    def send_request(self, data):
        to_send = pickle.dumps(data)
        for client in clients: 
            requests.post(f'http://{client["external_host"]}:{client["port"]}/api/pca/ld',
                data=to_send)

    def update(self,message):
        msg = pickle.loads(message)
        del message
        for key, val in msg.items():
           if key in self.sumLin:
               self.sumLin[key] += val[0]
               self.sumSq[key]  += val[1]
               self.cross[key]  += val[2]
           else:
               self.sumLin[key] = val[0]
               self.sumSq[key]  = val[1]
               self.cross[key]  = val[2]
        self.counter += 1
        if self.counter == self.num_clients: # Everyone reported
            message = {}
            for key, val in msg.items():
                corr_tot = corr(self.sumLin[key], self.sumSq[key],
                    self.cross[key])
                group  = store[key]
                tokeep = store["{}/PCA_passed".format(key)].value
                #positions = self.store["{}/PCA_positions".format(key)].value
                end = min(self.r1 + self.r0, len(tokeep))
                maf = store["{}/PCA_allele_freq".format(key)].value[
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
                                elif corr_tot[i,j] ** 2 > self.thresh:
                                    if maf[i] > (maf[j] * (1.0 + 
                                        Settings.kSmallEpsilon)):
                                        unfiltered[i] = False
                                    else:
                                        unfiltered[j] = False
                                    break
                    r2 = corr_tot[unfiltered,:][:,unfiltered]
                    if np.max(r2**2) < self.thresh:
                        break
                tokeep[self.r0:end][tokeep[self.r0:end]] = unfiltered
                pca_passed = store["{}/PCA_passed".format(key)]
                pca_passed[:] = tokeep
                end = min(self.r1 + self.r0 + self.r2, len(tokeep))
                if self.r2 + self.r0 >= len(tokeep) - 1:
                    message[key] = "END"
                else:
                    message[key] = tokeep[self.r0 + self.r2:end]
            self.send_request(message)
            self.counter = 0
            self.r0 += self.r2
            self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
#            self.tqdm.update(self.r2)

def report_pos(client_name):
    for chrom in store:
        if chrom == "meta":
            continue
        data = {}
        data[chrom] = store[f"{chrom}/PCA_passed"]
        msg = pickle.dumps(data)
        requests.post(f'http://{client["external_host"]}:{client_name["port"]}/api/pca/pcapos', msg)

