# stdlib
import logging
import pickle
import re
import os
import pdb
import time

# third party lib
import requests
import h5py
import numpy as np
from scipy.stats import chi2
chi2sf = chi2.sf

# internal lib
from lib.settings import Settings, Options, PCAFilterNames, Commands
from lib.client_registry import Registry
from lib.utils import write_or_replace


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


class LogisticAdmm(object):
    __instance = None
    def __init__(self, npcs, active, max_iters=20):
        if LogisticAdmm.__instance is not None:
            return
        else:
            self.estimates = {}
            self.iters = {}
            self.chroms = [v for v in store.keys() if v != 'meta']
            self.active_chrom = self.chroms[:active]
            self.finished = False
            self.likelihood = []
            self.max_iters = 5
            self.normalization_stats = None
            self.base_likelihood = 0
            self.nconnections = len(clients)
            LogisticAdmm.__instance = self
            self.send_request({}, "letTheDogsOut")

    def activate_chrom(self, chrom):
        msg = {"Estimated": "Small"}
        self.send_request(msg, "estimate")

    def make_chrom_active(self, chrom):
        self.chroms.remove(chrom)
        msg = {"Estimated":chrom}
        self.send_request(msg, "estimate")

    @staticmethod
    def get_instance(npcs, active):##TODO this is dangerous. If num_clients or win_size changes
        if LogisticAdmm.__instance is None:
            LogisticAdmm(npcs, active)
        return LogisticAdmm.__instance

    def update_stats(self, data):
        #TODO this should just be it's own object but for now it's faster to just add it here
        logging.info(f"Constructing feature's matrix!")
        data = pickle.loads(data)
        if self.normalization_stats is None:
            self.normalization_stats = data.copy()
            self.normalization_stats["iter"] = self.nconnections - 1
        else:
            self.normalization_stats["Sums"] += data["Sums"]
            self.normalization_stats["SS"] += data["SS"]
            self.normalization_stats["N"] += data["N"]
            self.normalization_stats["iter"] -= 1
            if self.normalization_stats["Indx"] != data["Indx"]:
                raise NameError("""Index of quantitative variables
                    does not match""")
            if self.normalization_stats["iter"] == 0:
                n = float(self.normalization_stats["N"])
                mu = self.normalization_stats["Sums"]/n
                pdb.set_trace()
                msg = {"Indx": self.normalization_stats["Indx"], "Means": mu,
                    "SD": np.sqrt(self.normalization_stats["SS"]/n - mu**2)}
                self.send_request(msg, "adjust")

                for chrom in self.active_chrom:
                    self.make_chrom_active(chrom)


    def send_request(self, data, subtask):
        to_send = pickle.dumps(data)
        for client in clients: 
            requests.post(f'http://{client["external_host"]}:{client["port"]}/api/asso/{subtask}',
                data=to_send)

    def update(self, message):
        message = pickle.loads(message)
        z_hat = message["VALS"]
        model = message["Estimated"]
        logging.info("Updating Estimate from {model}")
        self.estimates = self.update_estimate(z_hat, self.estimates)
    
    def association_finished(self):
        return self.finished
    
    def set_clients_state(self, state):
        for client in clients:
            Registry.get_instance().set_client_state(client['name'], state)

    def update_estimate(self, z_hat, model):
            if model in self.estimates:
                prev = self.estimates[model]
                if prev[1] == 1:
                    beta = (prev[0] + z_hat)/(self.nconnections)
                    self.iters[model] += 1
                    if self.iters[model] == self.max_iters:
                        write_or_replace(self.store, model + "/coefs", beta)
                        del self.estimates[model]
                        self.finished.add(model)
                        chroms = [key for key in self.store if key != 'meta']
                        if len(self.finished) == len(chroms) + 1:
                            self.set_clients_state("ASSO_DONE")
                            logging.info(f"Updating chr{chrom}. And now we are done!")
                            self.finished = True
                            self.initialize_pval_computation()
                    self.estimates[model] = [beta, self.nconnections - 1]
                    logging.info(f"Finished iteration{self.iters[model]} on chrom {model}")
                    msg = {"Estimated": model, "VALS": beta}
                    self.send_request(msg, "estimate")
                else:
                    self.estimates[model] = [prev[0] + z_hat , prev[1]-1]

            else: # Not in dictionary yet
                self.estimates[model] = [z_hat, self.nconnections - 1]
                self.iters[model] = 1

    def initialize_pval_computation(self):
        self.chroms = set([v for v in store.keys() if v != 'meta'])
        self.send_coef("Small", {})

    def send_coef(self, model, msg):
        coef = self.store[f"{model}/coefs"]
        msg["Estimated"] = model
        msg["Coef"] = coef
        self.send_request(msg, "coef")

    def update_pval(self, message):
        message = pickle.load(message)
        model = msg["Estimated"]
        if model in self.likelihood:
            prev  = self.likelihood[model]
            val = message["estimate"]
            if prev[1] == 1: 
                ell = prev[0] + val
                if model == "Small":
                    self.base_likelihood = ell
                    self.store.attrs["baselikelihood"] = ell
                else:
                    ell -= self.base_likelihood
                    pval = chi2sf(2*ell,1)
                    write_or_replace(self.store, f"{model}/ell", ell)
                    write_or_replace(self.store, f"{model}/llr_pval", pval)
                    del self.likelihood[model]
                new_model = self.chroms.pop()
                self.send_coef(new_model, {})
            else:
                self.likelihood[model] =[prev[0] + val, prev[1]-1]
        else:
            self.likelihood[model] =[val, self.nconnections-1]

