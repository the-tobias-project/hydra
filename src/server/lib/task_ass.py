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

