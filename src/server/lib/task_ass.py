# stdlib
import logging
import pickle
import os
import time

# third party lib
import h5py
import numpy as np
from scipy.stats import chi2
from numpy.linalg import _umath_linalg
from flask import current_app as app


# internal lib
from lib.settings import Settings, Options, PCAFilterNames, Commands
from lib.client_registry import Registry
from lib.utils import write_or_replace
import lib.networking as networking
from server.lib.plots import manhattan_plot
import pdb


#from  sklearn.linear_model import LogisticRegression as lr

storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()
chi2sf = chi2.sf
lstsq = _umath_linalg.lstsq_m
dot = np.dot


class LogisticAdmm(object):
    __instance = None
    def __init__(self, npcs, active, max_iters=50, algo="newton"):
        if LogisticAdmm.__instance is not None:
            return
        else:
            self.estimates = {}
            self.iters, self.accumulant = {}, {}
            self.chroms = [v for v in store.keys() if v != 'meta']
            self.active_chroms = ["Small"]
            self.finished = {}
            self.likelihood = {}
            self.max_iters = max_iters  * 5
            self.normalization_stats = None
            self.base_likelihood = 0
            self.nconnections = len(clients)
            LogisticAdmm.__instance = self
            self.send_request({}, "initialize")
            self.beta = 0
            self.algo=algo# "newton" for admm
            self.Hess, self.Gradients, self.Diags, self.Vals = {}, {}, {}, {}
            self.converged, self.linesearch_convergence = {}, {}
            self.threshold = .005
            self.linesearch_iter, self.fchanges = {}, {}
            self.scratch_likelihoods = {}
            self.alpha, self.BETA = .1, .5

            self.covars = {}

    def activate_chrom(self, chrom):
        msg = {"Estimated": "Small"}
        self.send_request(msg, "estimate")

    def make_chrom_active(self, chrom):
        self.finished[chrom] = False
        msg = {"Estimated":chrom}
        if self.algo != "admm":
            L = store[f"{chrom}/positions"].shape[0]
            starter = self.estimates["Small"]
            warm_start = np.empty((starter.shape[0]+1, 1))
            warm_start[1] = 0
            warm_start[0] = starter[0]
            warm_start[2:] = starter[1:]
            warm_start = np.tile(warm_start, (L,1,1))
            msg["VALS"] = warm_start
            self.estimates[chrom] = warm_start
        self.send_request(msg, "estimate")

    @staticmethod
    def get_instance(args, active):
        if LogisticAdmm.__instance is None:
            LogisticAdmm(args["ASSO_PCS"], active)
        return LogisticAdmm.__instance

    def update_stats(self, data):
        #TODO this should just be its own object but for now it's faster to just add it here
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
            #self.normalization_stats["data"] = np.vstack((self.normalization_stats["data"], data["data"]))
            #self.normalization_stats["ys"] = np.vstack((self.normalization_stats["ys"], data["ys"]))
            if self.normalization_stats["iter"] == 0:
                #cov = self.normalization_stats["data"]
                #Ys = self.normalization_stats["ys"]
                ind = self.normalization_stats["Indx"]

                n = float(self.normalization_stats["N"])

                mu = self.normalization_stats["Sums"]/n
               # pdb.set_trace()

               # cov[:,ind] -= mu
               # cov[:, ind] /= np.sqrt(self.normalization_stats["SS"]/n - mu**2)
               # ind = [i for i in range(cov.shape[1]) if i != 1]
               # model = lr(C=1e10, fit_intercept = False, tol=1e-7)
               # model.fit(cov[:,ind], Ys)
                logging.info(f"Total of {n} samples")
                msg = {"Indx": self.normalization_stats["Indx"], "Means": mu,
                    "SD": np.sqrt(self.normalization_stats["SS"]/n - mu**2)}
                self.send_request(msg, "adjust")

    def send_request(self, data, subtask):
        to_send = pickle.dumps(data)
        networking.message_clients(f"asso/{subtask}", data=to_send, env=app.config["ENV"])

    def update(self, message):
        message = pickle.loads(message)
        #if "H" in message:
        #    pdb.set_trace()
        z_hat = message["VALS"]
        model = message["Estimated"]
        logging.info(f"Updating Estimate from {model}")
        self.update_estimate(z_hat, model, message)

    def association_finished(self):
        return self.finished

    def set_clients_state(self, state):
        for client in clients:
            Registry.get_instance().set_client_state(client['name'], state)

    def update_estimate(self, z_hat, model, data):
        if model in self.estimates:
            if self.iters[model] >= self.max_iters: # this shouldn't happen but it does! WHy?
                logging.info(f"WHYYYYY {model}, {self.iters[model]}")
                return
            prev = self.estimates[model]
            #self.normalization_stats["data"] = np.vstack((self.normalization_stats["data"], data["cov"]))
            if prev[1] == 1:
                logging.info(f"Finished iteration {self.iters[model]} on chrom {model}")
                beta = (prev[0] + z_hat)/(self.nconnections)
                self.iters[model] += 1
                if model == "Small":
                    #print(np.linalg.norm(beta - self.beta))
                    print(np.sum(np.abs(beta - self.beta)))
                    self.beta = beta
                else:
                  print(np.linalg.norm(beta[7,:] - self.beta))
                  self.beta = beta[7,:]


                if self.iters[model] == self.max_iters:#TODO this shouldn't happen but why does it?
                    write_or_replace(store, f"meta/{model}/coef", beta)
                    if model != "Small":
                        del self.estimates[model]
                    self.active_chroms.remove(model)
                    if model == "Small":
                        self.max_iters = 15
                        self.estimates["Small"] = beta
                    #chroms = [key for key in store if key != 'meta']
                    #if len(self.finished) == len(chroms) + 1:
                    self.finished[model] = True
                    if not self.chroms:
                        self.set_clients_state("ASSO_DONE")
                        logging.info(f"We are done with association!")
                        self.initialize_pval_computation()
                    else:
                        chrom = self.chroms.pop()
                        self.active_chroms.append(chrom)
                        self.make_chrom_active(chrom)
                    return
                else:
                    self.estimates[model] = [np.zeros_like(beta), self.nconnections]
                    msg = {"Estimated": model, "VALS": beta, "Iter": self.iters[model]}
                    self.send_request(msg, "estimate")
            else:
                self.estimates[model] = [prev[0] + z_hat , prev[1]-1]

        else: # Not in dictionary yet
            #self.normalization_stats["data"] = data["cov"]

            self.estimates[model] = [z_hat, self.nconnections - 1]
            self.iters[model] = 1

    def newton_stats_update(self, msg, unconverged=None):
        data = pickle.loads(msg)
        model = data["Estimated"]
        if model in self.accumulant:
            #if self.iters[model] >= self.max_iters: # this shouldn't happen but it does! WHy?
            #    logging.info(f"WHYYYYY {model}, {self.iters[model]}")
            #    return
            self.Gradients[model] += data['g']
            self.Diags[model] += data['d']
            if model in self.iters:
                temp = self.Vals[model]
                temp[np.logical_not(self.converged[model])[:,0]] += data['v']
                self.Vals[model] = temp
            else:
                self.Vals[model] += data['v']
            self.Hess[model] += data['H']
            self.covars[model] = np.vstack((self.covars[model], data['covar']))
            self.accumulant[model] -= 1
        else:
            if model in self.Vals:
                temp = self.Vals[model]
                temp[np.logical_not(self.converged[model])[:,0]] = data['v']
                self.Vals[model] = temp
            else:
                self.Vals[model] = data['v']
            self.Gradients[model] = data['g']
            self.Diags[model] = data['d']
            self.Hess[model] = data['H']
            self.covars[model] = data['covar']
            self.accumulant[model] = self.nconnections - 1
        #if self.accumulant[model] == 0:
        #    pdb.set_trace()
        #    h = np.triu(self.Hess[model][25,:,:])
        #    h += h.T
        #    h += np.diag(self.Diags[model][51,:])
        finished_collecting = False
        if self.accumulant[model] == 0:
            finished_collecting = True
            del self.accumulant[model]
        return model, finished_collecting

    def newton_iter(self, model):
        TOL = 1e-7
        self.t = 1
        hs = self.Hess[model]
        gs = self.Gradients[model]
        ds = self.Diags[model]
        dfxs = np.zeros((gs.shape[0],1))
        af = store[f"{model}/allele_freq"].value
        x0 = self.estimates[model]
        if model not in self.converged:
            self.iters[model] = 1
            convergence_status = np.zeros((gs.shape[0], 1), dtype=bool)
            self.converged[model] = convergence_status.copy()
            old_not_converged = None
        else:
            self.iters[model] += 1
            old_not_converged = np.logical_not(self.converged[model])
            af = af[old_not_converged[:,0]]
            x0 = x0[old_not_converged[:,0]]
            convergence_status = np.zeros((x0.shape[0],1), dtype=bool)
        for i, g in enumerate(gs):
            if af[i] < self.threshold or 1-af[i] < self.threshold:
                convergence_status[i] = True
                continue
            h = np.diagflat(ds[i])
            if i % 2:
                hessian = np.triu(hs[i//2])
            else:
                hessian = np.tril(hs[i//2])
            hessian += hessian.T
            hessian += h
            dx = lstsq(-hessian, g[:,np.newaxis], 1e-5, signature='ddd->ddid')[0]
            dfx = dot(g.T, dx)
            if np.abs(dfx) < TOL:
                x0[i,:] += dx
                convergence_status[i] = True
                continue
            gs[i:i+1,:] = x0[i].T + self.t*dx.T
            dfxs[i,:] = dfx
            ds[i] = dx[:,0]
        if old_not_converged is not None:
            self.estimates[model][old_not_converged[:,0]] = x0
        else:
            self.estimates[model] = x0
        if np.prod(convergence_status) or self.iters[model] == self.max_iters:
            self.finished[model] = True
            write_or_replace(store, f"meta/{model}/newton_coef", self.estimates[model])
            af = store[f"{model}/allele_freq"].value
            arr = np.logical_not(np.logical_or(af < self.threshold, 1-af < self.threshold))
            est = self.estimates[model][arr]
            msg = {"Estimated": model, "conv": np.expand_dims(arr, axis=1), "x0": est[:,:,0]}
            self.send_request(msg, "query")
            del self.estimates[model]
            self.active_chroms.remove(model)
            if not self.chroms:
                self.set_clients_state("ASSO_DONE")
                logging.info("We are done with association")
                #COMPUTE PVALS
            else:
                chrom = self.chroms.pop()
                self.active_chroms.append(chrom)
                self.make_chrom_active(chrom)
            return
        else:
            if self.iters[model] > 1:
                temp = self.converged[model]
                temp[np.logical_not(temp)] = convergence_status[:,0]
                self.converged[model] = temp
            else:
                self.converged[model] = convergence_status
            gs = gs[np.logical_not(convergence_status)[:,0],:]
            msg = {"Estimated": model, "x0": gs, "conv": np.logical_not(self.converged[model])}
            self.send_request(msg, "query")
        self.linesearch_convergence[model] = convergence_status.copy()
        del self.Hess[model], self.Gradients[model]
        self.Diags[model] = ds[np.logical_not(convergence_status)[:,0],:]  # repurposed
        #self.fchanges[model] = dfxs[np.logical_not(convergence_status)]
        self.fchanges[model] = dfxs

    def collect_likelihoods(self, data):
        data = pickle.loads(data)
        model = data["estimated"]
        if model in self.scratch_likelihoods:
            self.scratch_likelihoods[model] += data['v']
            if self.linesearch_iter[model] == 1:
                if self.finished[model]:
                    write_or_replace(store, f"meta/{model}/newton_ell",
                        self.scratch_likelihoods[model])
                    write_or_replace(store, f"meta/{model}/newton_pval",
                        chi2sf(-2*self.scratch_likelihoods[model],1))
                    del self.scratch_likelihoods[model]
                    del self.Hess[model], self.Gradients[model]
                    del self.Diags[model], self.fchanges[model]
                    if all(value for value in self.finished.values()):
                        manhattan_plot(storePath, "manhattan_plot.png")
                else:
                    self.newton_test_new_point(model)
        else:
            self.scratch_likelihoods[model] = data['v']
            self.linesearch_iter[model] = self.nconnections
        self.linesearch_iter[model] -= 1

    def newton_test_new_point(self, model):
        vals = self.scratch_likelihoods[model].copy()
        del self.scratch_likelihoods[model]
        notConverged = np.logical_not(self.linesearch_convergence[model]).copy()
        notOverallConverged = np.logical_not(self.converged[model])
        old_vals = self.Vals[model][notOverallConverged]
        fs = old_vals + self.t * self.alpha * self.fchanges[model][notConverged]
        #fs = old_vals[notConverged[:,0]] + self.t * self.alpha * self.fchanges[model][notConverged[notOverallConverged]]
        #fs = self.Vals[model][notConverged] + self.t * self.alpha * self.fchanges[model]
        #[notOverallConverged[notConverged]]
        #fs = self.Vals[model][notConverged] + self.t * self.alpha * self.fchanges[model][notConverged[np.logical_not(self.converged[model])]]
        toUpdate = vals < fs[:, np.newaxis]
        nowConverged = notConverged.copy()
        nowConverged[notConverged] = np.logical_and(notOverallConverged[notOverallConverged], toUpdate[:,0])
        #accept = np.logical_and(toUpdate, notConverged np.logical_not(convergence_status[convergence_status][:,np.newaxis]))
        temp = self.estimates[model][np.logical_not(self.converged[model][:,0])]
        temp[toUpdate[:,0],:,0] += self.t * self.Diags[model][toUpdate[:,0]]
        self.estimates[model][np.logical_not(self.converged[model][:,0])] = temp
        #self.estimates[model][nowConverged[:,0],:,0] += self.t * self.Diags[model][toUpdate[:,0]]
        #self.estimates[model][np.logical_not(self.linesearch_convergence[model])[:,0],:][accept[:,0],:][:,:,0] += self.t * self.alpha * self.Diags[model][toUpdate[:,0]]
        # Update function value
        likelihoods = self.Vals[model][np.logical_not(self.converged[model][:,0])]
        likelihoods[toUpdate[:,0]] = vals[toUpdate[:,0]]
        self.Vals[model][np.logical_not(self.converged[model][:,0])] = likelihoods
        notConverged[notConverged] = np.logical_not(toUpdate[:,0])
        if not np.prod(toUpdate): # not all approved. Start another line search
            logging.info(f"chrom: {model} | t: {self.t} | percent converged line search: {(1-np.mean(notConverged))*100}%")
            self.t *= self.BETA
            x0s = self.estimates[model][notConverged[:,0]][:,:,0] + (
                self.t * self.Diags[model][np.logical_not(toUpdate)[:,0]])
            msg = {"Estimated": model, "x0": x0s, "conv": notConverged}
            self.send_request(msg, "query")
            temp = self.Diags[model]
            self.Diags[model] = temp[np.logical_not(toUpdate[:,0]), :]
            self.linesearch_convergence[model] = np.logical_not(notConverged)
        else: # All approved, start a new iteration of main algo
            msg = {"Estimated": model, "unconv": np.logical_not(self.converged[model]),
                'VALS': self.estimates[model][np.logical_not(self.converged[model])[:,0],:] }
            self.send_request(msg, "estimate")
            del self.Diags[model]

    def initialize_pval_computation(self):
        self.chroms = set([v for v in store.keys() if v != 'meta'])
        self.send_coef("Small", {})

    def send_coef(self, model, msg):
        coef = store[f"meta/{model}/coef"].value
        msg["Estimated"] = model
        msg["Coef"] = coef
        self.send_request(msg, "coef")

    def update_pval(self, message):
        def _share_new_ll():
            if len(self.chroms) >= 1:
               new_model = self.chroms.pop()
               self.send_coef(new_model, {})

        message = pickle.loads(message)
        model = message["Estimated"]
        if model == "Small":
            _share_new_ll()
            return
        val = message["estimate"]
        if model in self.likelihood:
            prev  = self.likelihood[model]
            if prev[1] == 1:
                #pdb.set_trace()
                ell = prev[0] + val
                pval = chi2sf(-2*ell,1)
                write_or_replace(store, f"meta/{model}/ell", ell)
                write_or_replace(store, f"meta/{model}/pval", pval)
                del self.likelihood[model]
                _share_new_ll()
            else:
                self.likelihood[model] =[prev[0] + val, prev[1]-1]
        else:
            self.likelihood[model] =[val, self.nconnections-1]

