# stdlib
import logging
import pickle
import os
import time

# third party lib
import h5py
import numpy as np
from flask import current_app as app
from scipy.sparse.linalg import eigsh as eig

# internal lib
from lib.settings import Settings, PCAFilterNames, Commands, ServerHTTP
from lib.client_registry import Registry
from lib import networking
from . import task_qc
from lib.corr import corr


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


def ready_to_decompose():
    chroms = [v for v in store if v != "meta"]
    if "meta" not in store:
        return False
    meta = store["meta"]
    print(chroms)
    for ch1 in chroms:
        for ch2 in chroms:
            if "{}_{}".format(ch1, ch2) not in meta and "{}_{}".format(ch2, ch1) not in meta:
                return False
    return True


def filtered():
    chrom = [key for key in store if key != "meta"][0]
    if "PCA_passed" in store[chrom] :
        return True
    return False


def start_pca_filters(filters):
    logging.info("Initiating PCA filters...")
    qc_subset_filters = {key: value for key, value in filters.items()
                         if key != PCAFilterNames.PCA_LD}
    if qc_subset_filters:  # Perform the qc like filters
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
            self.thresh = thresh
            CovarianceAggregator.__instance = self

    @staticmethod
    def get_instance(num_clients, win_size):##TODO this is dangerous. If num_clients or win_size changes
        if CovarianceAggregator.__instance is None:
            CovarianceAggregator(num_clients, win_size)
        return CovarianceAggregator.__instance

    def send_request(self, data, params=None):
        to_send = pickle.dumps(data)
        if params is None:
            networking.message_clients("pca/ld", env=app.config["ENV"], data=to_send)
        else:
            networking.message_clients("pca/ld", env=app.config["ENV"], data=to_send, args=params)

    def update(self, message):
        msg = pickle.loads(message)
        del message
        for key, val in msg.items():
            if key in self.sumLin:
                self.sumLin[key] += val[0]
                self.sumSq[key] += val[1]
                self.cross[key] += val[2]
            else:
                self.sumLin[key] = val[0]
                self.sumSq[key] = val[1]
                self.cross[key] = val[2]
        self.counter += 1
        if self.counter == self.num_clients: # Everyone reported
            message = {}
            for key, val in msg.items():
                corr_tot = corr(self.sumLin[key], self.sumSq[key],
                                self.cross[key])
                group = store[key]
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
                        if not snp1:  # already filtered
                            continue
                        else:
                            for j in range(i+1, n):
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

class Position_reporter(object):
    __instance = None

    def __init__(self):
        if Position_reporter.__instance is not None:
            return
        else:
            self.incrementor = len(clients)
            Position_reporter.__instance = self

    @staticmethod
    def get_instance():##TODO this is dangerous. If num_clients or win_size changes
        if Position_reporter.__instance is None:
            Position_reporter()
        return Position_reporter.__instance

    def report_pos(self):
        if self.incrementor > 1:
            self.incrementor -= 1
        else:
            for chrom in store:
                    if chrom == "meta":
                        continue
                    data = {}
                    data[chrom] = store[f"{chrom}/PCA_passed"].value
                    msg = pickle.dumps(data)
                    networking.message_clients("pca/pcapos", data=msg, env=app.config["ENV"])
            time.sleep(ServerHTTP.wait_time) #TODO fix this hack. Give time for clients to finish
            networking.message_clients("pca/cov", data=msg, env=app.config["ENV"])




def report_pos(client_names=None):
    for chrom in store:
        if chrom == "meta":
            continue
        data = {}
        data[chrom] = store[f"{chrom}/PCA_passed"].value
        msg = pickle.dumps(data)
        networking.message_clients("pca/pcapos", data=msg, env=app.config["ENV"])
    time.sleep(ServerHTTP.wait_time) #TODO fix this hack. Give time for clients to finish
    networking.message_clients("pca/cov", data=msg, env=app.config["ENV"])

def store_covariance(client_name, data):
    # store the chunks in the store. build on them
    msg = pickle.loads(data)
    ch1 = msg["CH1"]
    ch2 = msg["CH2"]
    logging.info("dealing with {}_{}".format(ch1, ch2))
    if "meta" not in store:
        store.create_group("meta")
    group = store["meta"]
    #if 'metadata' in msg:
    #    metadata = msg['metadata']
    #    if 'namespace' in metadata:
    #        namespace = metadata['namespace']
    #        logging.info(f'buildCov() namespace: {namespace}')
    # iter_number = msg['curr_iter']
    # logging.info(f'Reading iteration number {iter_number}')
    mat = msg["MAT"]
    cov_name = "{}_{}".format(ch1, ch2)
    logging.info(cov_name)
    if cov_name in group:
        #mat += group[cov_name].value
        stored = group[cov_name]
        stored[:,:] += mat
    else:
        group.create_dataset(cov_name, data=mat)
    logging.info(cov_name)
    if "E" in msg:
        logging.info("Finished storing covariances")
        eigenDecompose(n_components=10)


def eigenDecompose(n_components=10):
    chroms = sorted([int(v) for v in store.keys() if v != 'meta'])
    cov_size = 0
    meta = store["meta"]
    if "Vs" not in meta:
        for chrom in chroms:
            cov_size += meta["{}_{}".format(chrom, chrom)].shape[0]
        logging.info("Starting covariance matrix of size {} x {}".format(
            cov_size, cov_size))
        cov = np.empty((cov_size, cov_size), dtype=np.float32)
        i_old = 0
        for chrom1 in chroms:
            j_old = 0
            for chrom2 in chroms:
                if chrom2 > chrom1:
                    continue
                cov_name = "{}_{}".format(chrom1, chrom2)
                if cov_name in meta:
                    pcov = meta[cov_name].value
                    print(cov_name, pcov.shape)
                    cov[i_old:i_old+pcov.shape[0], j_old:j_old+pcov.shape[1]] = pcov
                    cov[j_old:j_old+pcov.shape[1], i_old:i_old+pcov.shape[0]] = pcov.T
                    j_old += pcov.shape[1]
            i_old += pcov.shape[0]
        cov /= (cov.shape[0])
        sigma, v = eig(cov, k=n_components, ncv=3*n_components)  # , maxiter=20*cov.shape[0])
        sigma, v = zip(*sorted(zip(sigma, v.T), reverse=True))
        v = np.array(v)
        sigma = np.array(sigma)
        sigma[sigma < 0] = 0
        print(sigma)
        meta.create_dataset('Sigmas', data = sigma)
        meta.create_dataset('Vs', data = v)
    else:
        logging.info("Eigenvalue decomposition is already done")
        sigma = meta["Sigmas"].value
        v = meta["Vs"].value
    sigma = np.sqrt(sigma) * np.sqrt(v.shape[0])
    inv_sigma = sigma.copy()
    inv_sigma[inv_sigma>0] = 1 / inv_sigma[inv_sigma > 0]
    print(chroms)
    msg = {"ISIG": inv_sigma, "V": v, "CHROMS": chroms}
    msg = pickle.dumps(msg)
    networking.message_clients("pca/eig", data=msg, env=app.config["ENV"])
