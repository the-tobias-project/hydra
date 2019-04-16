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

# internal lib
from lib.settings import Settings, Options, PCAFilterNames, Commands
from lib.client_registry import Registry
from lib.utils import write_or_replace


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()
chi2sf = chi2.sf


class LogisticAdmm(object):
    __instance = None
    def __init__(self, npcs, active, max_iters=3):
        if LogisticAdmm.__instance is not None:
            return
        else:
            self.estimates = {}
            self.iters = {}
            self.chroms = [v for v in store.keys() if v != 'meta']
            self.chroms = ["22"] #TODO change this
            self.active_chroms = ["Small"]
            self.finished = False
            self.likelihood = {}
            self.max_iters = max_iters
            self.normalization_stats = None
            self.base_likelihood = 0
            self.nconnections = len(clients)
            LogisticAdmm.__instance = self
            self.send_request({}, "initialize")

    def activate_chrom(self, chrom):
        msg = {"Estimated": "Small"}
        self.send_request(msg, "estimate")

    def make_chrom_active(self, chrom):
        msg = {"Estimated":chrom}
        self.send_request(msg, "estimate")

    @staticmethod
    def get_instance(npcs, active):
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
                msg = {"Indx": self.normalization_stats["Indx"], "Means": mu,
                    "SD": np.sqrt(self.normalization_stats["SS"]/n - mu**2)}
                self.send_request(msg, "adjust")

    def send_request(self, data, subtask):
        to_send = pickle.dumps(data)
        for client in clients:
            requests.post(f'http://{client["external_host"]}:{client["port"]}/api/asso/{subtask}',
                data=to_send)

    def update(self, message):
        message = pickle.loads(message)
        z_hat = message["VALS"]
        model = message["Estimated"]
        logging.info(f"Updating Estimate from {model}")
        self.update_estimate(z_hat, model)

    def association_finished(self):
        return self.finished

    def set_clients_state(self, state):
        for client in clients:
            Registry.get_instance().set_client_state(client['name'], state)

    def update_estimate(self, z_hat, model):
        if model in self.estimates:
            if self.iters[model] >= self.max_iters: # this shouldn't happen but it does! WHy?
                logging.info(f"WHYYYYY {model}, {self.iters[model]}")
                return
            prev = self.estimates[model]
            if prev[1] == 1:
                logging.info(f"Finished iteration {self.iters[model]} on chrom {model}")
                beta = (prev[0] + z_hat)/(self.nconnections)
                self.iters[model] += 1
                self.estimates[model] = [beta, self.nconnections]
                if self.iters[model] == self.max_iters:#TODO this shouldn't happen but why does it?
                    write_or_replace(store, f"meta/{model}/coef", beta)
                    del self.estimates[model]
                    self.active_chroms.remove(model)
                    #chroms = [key for key in store if key != 'meta']
                    #if len(self.finished) == len(chroms) + 1:
                    if not self.chroms:
                        self.set_clients_state("ASSO_DONE")
                        logging.info(f"We are done with association!")
                        self.finished = True
                        self.initialize_pval_computation()
                    else:
                        chrom = self.chroms.pop()
                        self.active_chroms.append(chrom)
                        self.make_chrom_active(chrom)
                else:
                    msg = {"Estimated": model, "VALS": beta}
                    self.send_request(msg, "estimate")
            else:
                self.estimates[model] = [prev[0] + z_hat , prev[1]-1]

        else: # Not in dictionary yet
            self.estimates[model] = [z_hat, self.nconnections - 1]
            self.iters[model] = 1

    def initialize_pval_computation(self):
        self.chroms = set([v for v in store.keys() if v != 'meta'])
        self.chroms = ["22"] #TODO change this
        self.send_coef("Small", {})

    def send_coef(self, model, msg):
        coef = store[f"meta/{model}/coef"].value
        msg["Estimated"] = model
        msg["Coef"] = coef
        pdb.set_trace()
        self.send_request(msg, "coef")

    def update_pval(self, message):
        message = pickle.loads(message)
        model = message["Estimated"]
        val = message["estimate"]
        if model in self.likelihood:
            prev  = self.likelihood[model]
            if prev[1] == 1: 
                ell = prev[0] + val
                if model == "Small":
                    self.base_likelihood = ell
                    store.attrs["baselikelihood"] = ell
                else:
                    ell -= self.base_likelihood
                    pval = chi2sf(2*ell,1)
                    write_or_replace(store, f"meta/{model}/ell", ell)
                    write_or_replace(store, f"meta/{model}/llr_pval", pval)
                    del self.likelihood[model]
                if len(self.chroms) > 1:
                    new_model = self.chroms.pop()
                    self.send_coef(new_model, {})
            else:
                self.likelihood[model] =[prev[0] + val, prev[1]-1]
        else:
            self.likelihood[model] =[val, self.nconnections-1]

