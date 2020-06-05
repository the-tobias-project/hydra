# stdlib
import logging
import os
import pickle
import time

# third party lib
import h5py
import numpy as np
from lib import networking
from flask import current_app as app

# internal lib
from lib.corr import hweP
from lib.settings import Settings, Commands
from lib.utils import write_or_replace
from lib.client_registry import Registry
from server.lib import plots
from lib import tasks


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
TIME = time.time()


class Echo(object):
    __instance = None

    def __init__(self, count):
        if Echo.__instance is not None:
            return
        else:
            Echo.count = count
            Echo.t = time.time()
            Echo.echos_left = count
            Echo.__instance = self
            networking.message_clients("echo", env=app.config["ENV"])
    
    @staticmethod
    def get_instance(count):
        if Echo.__instance is None: 
            Echo(count)
        return Echo.__instance

    def echo(self, client_name):
        instances = Registry.get_instance()
        instances.set_client_state(client_name, Commands.ECHO)
        if instances.num_clients_in_state(Commands.ECHO) == len(instances):
            Echo.echos_left -= 1
            for client in instances.list_clients():
                instances.set_client_state(client["name"], None)
        if not Echo.echos_left:
            Echo.__instance = None # remove the instance (essentially)
            avg_t = (time.time() - Echo.t)/ Echo.count
            networking.message_clients("End_echo", env=app.config["ENV"], data=avg_t)
            return avg_t
        else:
            networking.message_clients("echo", env=app.config["ENV"])
        return None


def start_init_task():
    global TIME
    TIME = time.time()
    tr = tasks.TaskReg.get_instance()
    tr.set_up_task(task=Commands.INIT, subtask="Start")
    for client in Registry.get_instance().list_clients():
        Registry.get_instance().set_client_state(client['name'], Commands.INIT)
    #networking.message_clients("init", env=app.config["ENV"])


def store_positions(data, client_name):
    data = pickle.loads(data)
    chrom = data['CHROM']
    positions = data['POS']
    dsetname = "{}/positions".format(chrom)
    write_or_replace(store, dsetname, positions, np.uint32)
    logging.info(f"{client_name} has {len(positions)} loci in chromosome {chrom}.")


def store_counts(data, client_name):
    message = pickle.loads(data)
    n = message["n"]
    if "START" in message:
        if "N" not in store.attrs:
            store.attrs["N"] = 0
            logging.info('Storing counts.')
        store.attrs["N"] += n
    chrom = message["CHROM"]
    size = len(message["COUNTS"])
    dsetname = "{}/counts".format(chrom)
    if dsetname not in store:
        dset = store.require_dataset(dsetname, (size, 4), dtype=np.int64)
    else:
        dset = store[dsetname]
    counts = message["COUNTS"]
    homo_ref = n - np.sum(counts, axis=1)[:, np.newaxis].astype(
        np.int64)
    dset[:] += np.hstack((homo_ref, counts))
    if "END" in message:
        Registry.get_instance().set_client_state(client_name, Commands.INIT_STATS)
        if Registry.get_instance().num_clients_in_state(Commands.INIT) == 0:
            logging.info('Transfering QC summary stats.')
            count_stats()


def count_stats():
    N = float(store.attrs["N"])
    task = "INIT"
    clients = Registry.get_instance().list_clients()
    for chrom in store.keys():
        counts_dset = store["{}/counts".format(chrom)].value
        missing_rate = counts_dset[:, 3] / float(N)
        store.create_dataset(f"{chrom}/missing_rates", data=missing_rate)
        af = (counts_dset[:, 2] * 2 + counts_dset[:, 1]).astype(float)
        af /= (np.sum(counts_dset[:, :3], axis=1) * 2).astype(float)
        # af = np.minimum(af, 1-af)
        store.create_dataset("{}/allele_freq".format(chrom), data=af)
        # var = counts_dset[:,0] * (2*af)**2
        # var += counts_dset[:,1] * (1-2*af)**2
        # var += counts_dset[:,2] * (2-2*af)**2
        # var /= (N-counts_dset[:,3]) # 2*af*(1-af)
        var = 2 * af * (1 - af)
        store.create_dataset("{}/var".format(chrom), data=var)
        hwe = hweP(counts_dset[:, :3].astype(np.int32), 1, 0)
        # Need to Recompile HWEP with uint32
        msg = {
            "TASK": task,
            "SUBTASK": "STATS",
            "CHROM": chrom,
            "HWE": hwe.tolist(),
            "MISS": missing_rate.tolist(),
            "AF": af.tolist(),
            "VAR": var.tolist()
        }
        #msg = pickle.dumps(msg)
        tr = tasks.TaskReg.get_instance()
        tr.set_up_task(task=Commands.INIT, subtask="stats", other={"data": msg})
        #networking.message_clients("init/stats", env=app.config["ENV"], data=msg)
        store.create_dataset("{}/hwe".format(chrom), data=hwe)
    for client in clients:
        Registry.get_instance().set_client_state(client['name'], "DONE_INIT")

    logging.info("Done with initialization")
    make_plots("QC_pre_filter.png")
    logging.info(f"Initialization took roughly {time.time()-TIME:.1f} seconds.")


def make_plots(outname):
    logging.info("Generating QC plots")
    plots.qc_plots(storePath, outname)
