# stdlib
import os
import pickle
import sys
import time

# third party lib
import h5py
import numpy as np
from plinkio import plinkfile
from sklearn.metrics import log_loss

# internal lib
from client.lib import shared
from lib import HTTPResponse
from lib.utils import write_or_replace
from lib.optimizationAux import *
from lib.opt import minimize_lbfgs

#TODO clean up Us etc 

class LogisticAdmm(object):
    __instance = None
    def __init__(self, covar_numbers, npcs, client_config, threshold=.005):
        if LogisticAdmm.__instance is not None:
            return 
        else:
            pfile = client_config["plinkfile"]
            store_name = shared.get_plink_store(pfile)
            self.store = h5py.File(store_name, 'a')
            self.threshold = threshold
            self.client_config = client_config
            self.include_mask = []
            self.load_Y(covar_numbers[-1], pfile)
            self.load_covar_mat(npcs, covar_numbers[:-1], pfile)
            self.warm_start = 0
            self.previous_estimates = {}
            self.prev_cov_estimate = None
            self.previous_Us = {}
            self.flipped_covar = None
            self.baseline_likelihood = {}
            LogisticAdmm.__instance = self

    @staticmethod
    def get_instance(covar_numbers, npcs, client_config):
        if LogisticAdmm.__instance is None:
            LogisticAdmm(covar_numbers, npcs, client_config)
        return LogisticAdmm.__instance
    
    def load_Y(self, index, fname):
        name = shared.get_covar_file(fname)
        y = np.loadtxt(name, usecols=[index], delimiter='\t', skiprows=1)
        n = y.shape[0]
        self.include_mask = np.array(y != -9)
        print (f"{np.sum(self.include_mask)} out of {n} individuals have the phenotype.")
        y = y[self.include_mask]
        n = int(np.sum(self.include_mask))
        self.Ys = np.sign(y - .5).reshape(n, 1)

    def load_covar_mat(self, npcs, other_covars_loc, fname):
        store = self.store
        #n = store.attrs["n"]
        include_mask = self.include_mask
        n = np.sum(include_mask)
        p = npcs + len(other_covars_loc) + 2 # intercept + snp
        #if self.threshold is None: #TODO determine this based on entire sample size
        #    self.threshold = 3*float(p)/n   # Don't run regression unless you have at least this AF
        self.covariates = np.empty((n, p))
        # Load up scores 
        self.covariates[:, 2:npcs+2] = store["meta/pca_u"].value[include_mask,:npcs]
        covar_file_name = shared.get_covar_file(fname)
        other_covars = np.loadtxt(covar_file_name, delimiter="\t", usecols=other_covars_loc,
            skiprows=1)
        self.covariates[:, npcs+2:] = other_covars[include_mask, :]
        self.covariates[:, 0] = 1
        other_covars = self.send_summary_to_standardize()
        
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
        msg = {"Indx": quant_covars, "Sums": sums, "SS": sumsq, "N": self.covariates.shape[0]
            #,"data": self.covariates, "ys": self.Ys  #TODO take this out. only for debugging
            }
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
        #temp = self.covariates.copy()
        #self.covariates = temp[self.include_mask, :].copy()
        self.covariates *= -self.Ys
        self.flipped_covar = True
        selfmessage = pickle.dumps({"Estimated": "Small"})
        self.update(selfmessage, client_config)

    def update(self, data, client_config):
        data = pickle.loads(data)
        self.client_config = client_config
        y = self.Ys
        model = data["Estimated"]
        print(f"Updating {model}")
        beta = None
        if "VALS" in data:
            beta = data["VALS"]
        if model == "Small":
            self.run_covar_regression(warm_start=beta)
        else:
            #self.run_logistic_regression(y, model, warm_start=beta)
            in_mask = None
            if "unconv" in data:
                in_mask = data["unconv"]
            if model not in self.baseline_likelihood:
                exclude = [1]
                indx = np.array([i for i in range(self.covariates.shape[1]) if i not in exclude])
                self.baseline_likelihood[model] = self.evaluate_estimate(
                    model, in_mask, beta[:,indx,0], exclude=exclude)
            self.run_newton_lr(y, model, warm_start=beta, unconverged=in_mask)

    def run_covar_regression(self, warm_start = None, rho=250.0, alpha=1.0):
        # instead of making many function calls, I'll separate the covariate and chromosome regressions
        model = "Small"
        ncov = self.covariates.shape[1]
        estimates = np.zeros((ncov-1, 1))
        idx = [i for i in range(ncov) if i != 1]
        covariates = self.covariates[:,idx]
        if self.prev_cov_estimate is not None:
            z_hat = self.prev_cov_estimate
            all_Us = self.previous_Us[model] + z_hat - warm_start
        else:
            all_Us = 0
        #print(covariates)
        if warm_start is None:
            estimates[:, 0] = other_newton(covariates, np.zeros((ncov-1)),
                np.zeros((ncov-1,)), rho, estimates[:,0], ncov-1)
            z_hat = estimates
        else:
            estimates[:,0] = other_newton(covariates,all_Us[:,0],
                warm_start[:,0], rho, z_hat[:, 0], ncov-1)
            z_hat = alpha * estimates + (1-alpha) * warm_start
        self.prev_cov_estimate = estimates
        self.previous_Us[model] = all_Us
        est = z_hat+all_Us
        msg = pickle.dumps({"VALS": est, "Estimated":"Small"})
        HTTPResponse.respond_to_server('api/tasks/ASSO/estimate', 'POST', msg,
            self.client_config['name'])

    def run_newton_lr(self, y, chrom=None, warm_start=None, unconverged=None):
        store = self.store
        print("starting with newtons")
        include_mask = self.include_mask
        n = np.sum(include_mask)
        y = y.reshape(n)
        covariates = self.covariates.copy()
        group = store[chrom]
        af = group["MAF"]
        positions = group["positions"]
        ncov = self.covariates.shape[1] 
        baselikelihood = self.baseline_likelihood[chrom]
        if unconverged is None: 
            L = len(positions)
        else:
            L = np.sum(unconverged)
            positions = positions[unconverged[:,0]]
            af = af[unconverged[:,0]]
            baselikelihood = baselikelihood[unconverged[:,0]]

        hessians = np.zeros((int(np.ceil(L/2)), ncov, ncov))
        diagonals = np.zeros((L, ncov))
        gradients = np.zeros((L, ncov))
        vals = np.zeros((L, 1))
        count = 0
        t = time.time()
        for i, position in enumerate(positions):
          #  if i == 10: #TODO
          #      totalT = (time.time()-t)/500
          #      print(f"Average time is {totalT} for n={n}, {t2/500}")
          #      print(f"Count is {count}")
          #      break
            if not i % 1000:
                print(f"{time.time()-t}")
                t = time.time()
                print (float(i/L))
            if af[i] < self.threshold or (1-af[i]) < self.threshold:
                continue
            val = group[str(position)].value[include_mask]
            # Dumb imputation. Hopefully your data is already imputed and this doesn't happen
            val[np.isnan(val)] = 0
            #ind = ~np.isnan(val)
            #covariates[ind,1] = val[ind] * -y[ind]
            covariates[:,1] = val * -y
            count+=1
            #submat = covariates[ind,:]
            h,diagonals[i], gradients[i], vals[i,0] = ltri_Hessians(
                #submat, warm_start[i,:,0], ncov, submat.shape[0], 0) #rho set to zero
                covariates, warm_start[i,:,0], ncov, covariates.shape[0], 0) #rho set to zero
            if i % 2:
                hessians[i//2,:,:] += h.T
            else:
                hessians[i//2,:,:] += h
        vals -= baselikelihood
        msg = pickle.dumps({"Estimated": chrom, "H": hessians, 'g':gradients,
          'd': diagonals, 'v': vals, "covar": covariates})
        HTTPResponse.respond_to_server('api/tasks/ASSO/hessians', 'POST', msg,
            self.client_config['name'])

    def cost(self, data):
        msg = pickle.loads(data)
        chrom = msg["Estimated"]
        mask = msg["conv"]
        x0 = msg["x0"]
        estimates = self.evaluate_estimate(chrom, mask, x0)
        estimates -= self.baseline_likelihood[chrom][mask[:,0]]
        msg = pickle.dumps({'estimated': chrom, 'v': estimates})
        HTTPResponse.respond_to_server('api/tasks/ASSO/valback', 'POST', msg,
            self.client_config['name'])

    def evaluate_estimate(self, chrom, mask, x0, exclude=None):
        store = self.store
        include_mask = self.include_mask
        L = x0.shape[0]
        group = store[chrom]
        positions = group["positions"].value
        if mask is not None:
            positions = positions[mask[:,0]]
        covariates = self.covariates.copy()
        if exclude is not None:
            include = np.array([k for k in range(covariates.shape[1]) if k not in exclude])
            covariates = covariates[:, include]
        n = int(np.sum(self.include_mask))
        y = self.Ys.reshape(n)
        vals = np.zeros((L,1))
        for i, position in enumerate(positions):
            val = group[str(position)].value[include_mask]
            # the dumb imputation reappears
            val[np.isnan(val)] = 0
            #ind = ~np.isnan(val)
            if exclude is None or 1 not in exclude: # Not excluding the genotype or intercept
                #covariates[ind, 1] = val[ind] * -y[ind]
                covariates[:, 1] = val * -y
            vals[i, 0] = function_values(covariates[:,:], x0[i,:].T)
        return vals

#    def ref_impute(self):
#        store = self.store
#        chroms = [k for k in store if k!="meta"]
#        for chrom in chroms:
#            group = store[chrom]
#            positions = group[f"{positions}"].value
#            for position in positions:
#                dset = group[position]
#                temp = dset.value
#                temp[np.isnan(temp)] = 0
#                dset[...] = temp

    def run_logistic_regression(self, y, chrom=None, warm_start=None, rho=250.0, alpha=1.00):
        store = self.store
        include_mask = self.include_mask
        n = int(np.sum(self.include_mask))
        y = y.reshape(n)
        covariates = self.covariates.copy()
        group = store[chrom]
        positions = group["positions"]
        ncov = self.covariates.shape[1] 
        estimates = np.zeros((len(positions), ncov))
        if warm_start is None:
            est = np.zeros(ncov)
            covar_estimates = self.prev_cov_estimate
            # boundary condition for the loop
            est[1] = 0
            est[0] = covar_estimates[0]
            est[2:] = covar_estimates[1:].ravel()
        af = group["MAF"]
        if chrom in self.previous_estimates:
            z_hat = self.previous_estimates[chrom]
            all_Us = self.previous_Us[chrom] + z_hat - warm_start
        else:
            all_Us = np.zeros((ncov))
            all_Us[0] = self.previous_Us["Small"][0]
            all_Us[2:] = self.previous_Us["Small"][1:].ravel()
        count = 0
        for i, position in enumerate(positions):
          #  if i == 10: #TODO
          #      totalT = (time.time()-t)/500
          #      print(f"Average time is {totalT} for n={n}, {t2/500}")
          #      print(f"Count is {count}")
          #      break
            if not i % 100:
                print(f"{time.time()-t}")
                t = time.time()
                t2 = 0
                print (float(i/len(positions)))
            if af[i] < self.threshold or (1-af[i]) < self.threshold:
                estimates[i,:] = np.nan
                continue
            else:
                val = group[str(position)].value[include_mask]
                ind = ~np.isnan(val)
                covariates[ind,1] = val[ind] * -y[ind]
                count+=1
                t3 = time.time()
                if warm_start is None:
                    #est = np.ascontiguousarray(estimates[:,i])
                    #estimates[:,i] = minimize_lbfgs(covariates[ind,:], np.zeros((ncov)), np.zeros((ncov,)), rho, est, ncov)
                    #estimates[i,:] = minimize_lbfgs(covariates[ind,:], np.zeros((ncov)), est, rho, est, ncov)
                    estimates[i,:] = other_newton(covariates[ind,:], all_Us, est, rho, est, ncov)
                    #est *= .2
                    #est += .8*estimates[i,:]
                    #est[1] = 0
#                    print(f"{estimates[i,:]} from myc")
#                    estimates[i,:] = bfgs_more_gutted(covariates[ind,:], np.zeros((ncov)), np.zeros((ncov,)), rho, est, ncov)
#                    print(f"{estimates[i,:]} from sci")
                    z_hat = alpha * estimates + (1-alpha) * est
                else:
                    #estimates[i,:] = bfgs_more_gutted(covariates[ind,:], all_Us[i,:],
                    #    warm_start[i,:], rho, z_hat[i, :], ncov)
                    #estimates[i,:] = minimize_lbfgs(covariates[ind,:], all_Us[i,:],
                    #    warm_start[i,:], rho, z_hat[i,:], ncov)
                    estimates[i,:] = other_newton(covariates[ind,:], all_Us[i,:],
                        warm_start[i,:], rho, z_hat[i,:], ncov)

                    z_hat = alpha * estimates + (1-alpha) * warm_start
              #  t2 += time.time()-t3
        self.previous_estimates[chrom] = estimates
        self.previous_Us[chrom] = all_Us 
        msg = pickle.dumps({"Estimated": chrom, "VALS": z_hat + all_Us})#, 'cov': covariates})
        HTTPResponse.respond_to_server('api/tasks/ASSO/estimate', 'POST', msg,
            self.client_config['name'])


    def send_likelihood(self, message):
        #TODO Important, if we are excluding missing values, we should recompute baseline every time
        message = pickle.loads(message)
        include_mask = self.include_mask
        model = message["Estimated"]
        coef = message["Coef"]
        coef = coef.T
        covariates = self.covariates
        #n = self.store.attrs['n']
        n = int(np.sum(self.include_mask))
        y = self.Ys.copy()
        ell = None
        if self.flipped_covar:
            self.covariates *= -y
            self.flipped_covar = False
        y += 1
        y /= 2
        if model == "Small":
            indx = [i for i in range(covariates.shape[1]) if i != 1]
            y_model = 1.0 / (1 + np.exp(-covariates[:,indx].dot(coef.T)))
            self.base_y_pred = y_model
            #ell = log_loss((y+1)/2, y_model, normalize=False, labels=[0,1])
        else:
            group = self.store[model]
            af = group["MAF"].value
            tokeep = np.logical_and(af>self.threshold, 1-af>self.threshold)
            positions = group["positions"].value
            ell = np.zeros((1,positions.shape[0]))
            for i, position in enumerate(positions):
                if not tokeep[i]:
                    ell[0,i] = np.nan
                else:
                    val = group[str(position)].value[include_mask]
                    ind = ~ np.isnan(val) #TODO impute or something? 
                    covariates[:,1] = val
                    y_model = 1.0 / (1 + np.exp(-covariates[ind, :].dot(coef[:,i].T)))
                    ell[0,i] = log_loss(y[ind], y_model, normalize=False, labels=[0,1])
                    ell[0,i] -= log_loss(y[ind], self.base_y_pred[ind], normalize=False, labels=[0,1])

        msg = pickle.dumps({"Estimated": model, "estimate": ell})
        HTTPResponse.respond_to_server('api/tasks/ASSO/pval', 'POST', msg,
            self.client_config['name'])
