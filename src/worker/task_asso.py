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
from optimizationAux import *

class LogisticAdmm(object):
    __instance = None
    def __init__(self, covar_numbers, npcs, threshold=None):
        if LdReporter.__instance is not None:
            return 
        else:
            LdReporter.__instance = self
            pfile = client_config["plinkfile"]
            store_name = shared.get_plink_store(pfile)
            self.store = h5py.File(store_name, 'a')
            self.load_Y(covar_numbers[-1])
            self.load_covar_mat(npcs, covar_numbers[:-1])
            self.threshold = None
            self.warm_start = 0
            self.previous_estimates = {}
            self.prev_cov_estimate = None
            self.client_config = None
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
        if self.threshold is None:
            self.threshold = 3*float(p)/n   # Don't run regression unless you have at least this AF
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
        self.client_config = client_config
        y = self.Ys
        model = data["Estimated"]
        if model == "Small":
            self.covariate_regression(y)
        else:
            self.run_logistic_regression(y, model)
            pass # perform snp regression

    def run_covar_regression(self, y, warm_start = None, rho=1.0, alpha=1.0):
        # instead of making many function calls, I'll separate the covariate and chromosome regressions
        store = self.store
        n = store.attrs["n"]
        y = y.reshape(n)
        ncov = self.covariates.shape[1]
        estimates = np.zeros((ncov-2, 1))
        estimates[:, 0] = 0
        idx = [i for i in range(self.covariates.shape[1]) if i != 1]
        covariates = self.covariates[:,idx]
        if self.prev_cov_estimate is not None:
            z_hat = self.prev_cov_estimate
            all_Us = self.prev_Us + z_hat - warm_start
        else:
            all_Us = 0
        if warm_start is None:
            estimates[:, 0] = bfgs_more_gutted(covariates, np.zeros((ncov)),
                np.zeros((ncov,)), rho, estimates[:,i], ncov)
            z_hat = alpha * estimates
        else:
            estimates[:,0] = bfgs_more_gutted(covariates,all_Us[:,i],
                warm_start[:,i], rho, z_hat[:, i], ncov)
            z_hat = alpha * estimates + (1-alpha) * warm_start
        self.prev_cov_estimate = estimates
        self.previous_Us = all_Us
        msg = pickle.dumps({"VALS": z_hat, "Estimated":"Small"})
        HTTPResponse.respond_to_server('api/tasks/ASSO/estimate', 'POST', msg,
            self.client_config['name'])

    def run_logistic_regression(self, y, chrom=None, warm_start=None, rho=1.0, alpha=1.0):
        store = self.store
        n = store.attrs["n"]
        y = y.reshape(n)
        covariates = self.covariates
        group = store[chrom]
        positions = group["positions"]
        ncov = self.covariates.shape[1]
        estimates = np.zeros((ncov, len(positions)))
        estimates[:, -1] = 0
        af = group["MAF"]
        if chrom in self.previous_estimates:
            z_hat = self.previous_estimates[chrom]
            all_Us = self.previous_Us[chrom] + z_hat - warm_start
        else:
            all_Us = 0
        for i, position in enumerate(positions):
            if af[i] < self.threshold or (1-af[i]) < self.threshold:
                estimates[:, i] = np.nan
                continue
            covariates[:,1] = group[str(position)].value * -y
            if warm_start is None:
                estimates[:,i] = bfgs_more_gutted(covariates, np.zeros((ncov)), np.zeros((ncov,)), rho, estimates[:,i], ncov)
                z_hat = alpha * estimates
            else:
                estimates[:,i] = bfgs_more_gutted(covariates, all_Us[:,i],
                    warm_start[:,i], rho, z_hat[:, i], ncov)
                z_hat = alpha * estimates + (1-alpha) * warm_start
        self.previous_estimates[chrom] = estimates
        self.previous_Us[chrom] = all_Us
        msg = pickle.dumps({"Estimated": chrom, "VALS": z_hat})
        HTTPResponse.respond_to_server('api/tasks/ASSO/estimate', 'POST', msg,
            self.client_config['name'])
