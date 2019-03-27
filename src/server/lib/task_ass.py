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

# internal lib
from lib.settings import Settings, Options, PCAFilterNames, Commands
from lib.client_registry import Registry


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


class LogisticAdmm(object):
    __instance = None
    def __init__(self, npcs, active, max_iters=20):
        if LogisticAdmm.__instance is not None:
            return
        else:
            self.num_clients = num_clients
            self.estimates = {}
            self.small_model = {}
            self.iters = {}
            self.chroms = [v for v in store.keys() if v != 'meta']
            self.activ_chrom = self.chroms[:active]
            LogisticAdmm.__instance = self

    def activate_chrom(self, chrom):
        msg = {}
        self.send_request(self,)

    def make_chrom_active(self, chrom):
        self.chroms.remove(chrom)
        msg = {"CHROM":chrom}
        self.send_request(msg)

    @staticmethod
    def get_instance(num_clients, npcs):##TODO this is dangerous. If num_clients or win_size changes
        if LogisticAdmm.__instance is None:
            LogisticAdmm(npcs, active, num_clients, win_size)
        return LogisticAdmm.__instance

    def send_request(self, data):
        to_send = pickle.dumps(data)
        for client in clients: 
            requests.post(f'http://{client["external_host"]}:{client["port"]}/api/asso/estimate',
                data=to_send)

    def update(self, message):
        message = pickle.loads(message)
        z_hat = message["VALS"]
        model = message["Estimated"]

            self.small_model = update_estimate(z_hat, self.small_model)
        else:
            self.estimates = self.update_estimate(z_hat, self.estimates)


    def update_estimate(self, z_hat, model):
            if model in self.estimates:
                prev = self.estimates[model]
                if prev[1] == 1:
                    beta = (prev[0] + z_hat)/(self.connections)
                    self.iters[model] += 1
                    if self.iters[model] == self.max_iters:
                        write_or_replace(self.store, model + "/coefs", beta)
                        del self.estimates[model]
                        self.finished.add(model)
                        chroms = [key for key in self.store if key != 'meta']
                        if len(self.finished) == len(chroms) + 1:
                            logging.info(f"Updating chr{chrom}. And now we are done!")
                    self.estimates[model] = [beta, self.connections - 1]
                    logging.info(f"Finished iteration{self.iters[model]} on chrom {model}")
                    msg = {"Estimated": model, "VALS": beta}
                    self.send_request(msg)
                else:
                    self.estimates[model] = [prev[0] + z_hat , prev[1]-1]

            else: # Not in dictionary yet
                self.estimates[model] = [z_hat, self.connections - 1]
                self.iters[model] = 1

