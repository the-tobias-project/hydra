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
from sklearn.metrics import log_loss

# internal lib
from client.lib import shared
from lib import HTTPResponse
from lib.utils import write_or_replace
from lib.optimizationAux import *

class LogisticAdmm(object):
    __instance = None
    def __init__(self, covar_numbers, npcs, client_config, threshold=None):
        if LogisticAdmm.__instance is not None:
            return 
        else:
            pfile = client_config["plinkfile"]
            store_name = shared.get_plink_store(pfile)
            self.store = h5py.File(store_name, 'a')
            self.threshold = threshold
            self.client_config = client_config
            self.load_Y(covar_numbers[-1], pfile)
            self.load_covar_mat(npcs, covar_numbers[:-1], pfile)
            self.warm_start = 0
            self.previous_estimates = {}
            self.prev_cov_estimate = None
            self.flipped_covar = None
            LogisticAdmm.__instance = self

    @staticmethod
    def get_instance(covar_numbers, npcs, client_config):
        if LogisticAdmm.__instance is None:
            LogisticAdmm(covar_numbers, npcs, client_config)
        return LogisticAdmm.__instance
    
    def load_Y(self, index, fname):
        name = shared.get_covar_file(fname)
        y = np.loadtxt(name, usecols=[index], delimiter='\t')
        n = y.shape[0]
        self.Ys = np.sign(y - .5).reshape(n, 1)

    def load_covar_mat(self, npcs, other_covars_loc, fname):
        store = self.store
        n = store.attrs["n"]
        p = npcs + len(other_covars_loc) + 2 # intercept + snp
        if self.threshold is None:
            self.threshold = 3*float(p)/n   # Don't run regression unless you have at least this AF
        self.covariates = np.empty((n, p))
        # Load up scores 
        scores = store["meta/pca_u"].value[:,:npcs]
        scores /= np.std(scores, axis = 0)  # This needs to be the global std
        self.covariates[:, 2:npcs+2] = scores  # already centered
        del scores
        covar_file_name = shared.get_covar_file(fname)
        other_covars = np.loadtxt(covar_file_name, delimiter="\t", usecols=other_covars_loc)
        self.covariates[:, npcs+2:] = other_covars
        other_covars = self.send_summary_to_standardize()
        self.covariates[:, 0] = 1
        #other_covars /= np.std(other_covars, axis = 0)  # This needs to be the global std
        #other_covars -= np.mean(other_covars, axis = 0)  # This needs to be the global std
        #self.covariates[:, 0] = -self.Ys.reshape(-1)
        #self.flipped_covar = True
    
    def send_summary_to_standardize(self):
        # This is a rough sketch, the assumption here is that if you are quantitative factor, then you 
        # exhibit more than 2 values within each silo. I never actually verify this with the current version
        # So that might be a good TODO for future implementations. 
        quant_covars = [ i for i in range(2,self.covariates.shape[1]) if 
            len(np.unique(self.covariates[:,i])) > 2 ]
        sums = np.sum(self.covariates[:,quant_covars], axis = 0)
        sumsq = np.sum(self.covariates[:,quant_covars]**2, axis=0)
        msg = {"Indx": quant_covars, "Sums": sums, "SS": sumsq, "N": self.covariates.shape[0]}
        msg = pickle.dumps(msg)
        HTTPResponse.respond_to_server('api/tasks/ASSO/adjust', 'POST', msg,
            self.client_config['name'])

    def global_standardize(self, data, client_config):
        data = pickle.loads(data)
        indx = data["Indx"]
        mu = data["Means"]
        sd = data["SD"]
        self.covariates[:, indx] -= mu
        self.covariates[:, indx] /= sd
        self.covariates *= -self.Ys
        self.flipped_covar = True
        selfmessage = pickle.dumps({"Estimated": "Small"})
        self.update(selfmessage, client_config)

    def update(self, data, client_config):
        data = pickle.loads(data)
        self.client_config = client_config
        y = self.Ys
        model = data["Estimated"]
        if model == "Small":
            self.run_covar_regression(y)
        else:
            self.run_logistic_regression(y, model)

    def run_covar_regression(self, y, warm_start = None, rho=1.0, alpha=1.0):
        # instead of making many function calls, I'll separate the covariate and chromosome regressions
        store = self.store
        n = store.attrs["n"]
        y = y.reshape(n)
        ncov = self.covariates.shape[1]
        estimates = np.zeros((ncov-1, 1))
        estimates[:, 0] = 0
        idx = [i for i in range(self.covariates.shape[1]) if i != 1]
        covariates = self.covariates[:,idx]
        if self.prev_cov_estimate is not None:
            z_hat = self.prev_cov_estimate
            all_Us = self.prev_Us + z_hat - warm_start
        else:
            all_Us = 0
        print(covariates)
        if warm_start is None:
            estimates[:, 0] = bfgs_more_gutted(covariates, np.zeros((ncov-1)),
                np.zeros((ncov-1,)), rho, estimates[:,0], ncov-1)
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

    def send_likelihood(self, message):
        message = pickle.loads(message)
        model = message["Estimated"]
        coef = message["Coef"]
        covariates = self.covariates
        y = self.Ys
        if self.flipped_covar:
            self.covariates *= -y.reshape(-1)
        if model == "Small":
            indx = [i for i in range(covariates.shape[1]) if i != 1]
            y_model = 1.0 / (1 + np.exp(-covariates[:,indx].dot(coef)))
            ell = log_loss((y+1)/2, y_model, normalize=False, labels=[0,1])
        else:
            group = self.store[chrom]
            positions = group["positions"]
            ell = np.empty(1,positions.shape[1])
            for i, position in enumerate(positions):
                covariates[:,1] = group[str(position)].value 
                y_model = 1.0 / (1 + np.exp(-covariates.dot(coef)))
                ell = log_loss((y+1)/2, y_model, normalize=False, labels=[0,1])
        msg = pickle.dumps({"Estimated": model, "estimate": ell})
        HTTPResponse.respond_to_server('api/tasks/ASSO/pval', 'POST', msg,
            self.client_config['name'])
