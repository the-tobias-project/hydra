# stdlib
import logging
import os
import pickle

# third party lib
import h5py
import numpy as np
import requests

# internal lib
from lib.corr import corr, hweP
from lib.settings import Settings, Commands
from lib.utils import write_or_replace
from lib.client_registry import Registry


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


def start_init_task():
    for client in clients:
        Registry.get_instance().set_client_state(client['name'], Commands.INIT)
        requests.post(f'http://{client["external_host"]}:{client["port"]}/api/init')


def store_positions(data):
    data = pickle.loads(data)
    chrom = data['CHROM']
    positions = data['POS']
    dsetname = "{}/positions".format(chrom)
    write_or_replace(store, dsetname, positions, np.uint32)
    logging.info("{} loci in chromosome {}.".format(len(positions), chrom))


def store_counts(data, client_name):
    message = pickle.loads(data)

    n = message["n"]
    logging.info('storing counts')
    if "START" in message:
        if "N" not in store.attrs:
            store.attrs["N"] = 0
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
        Registry.get_instance().set_client_state(client_name, 'INIT_DONE')
        if Registry.get_instance().num_clients_in_state(Commands.INIT) == 0:
            logging.info('Done getting init reports from clients')
            logging.info('Telling clients to store stats')
            count_stats()


def count_stats():
    N = float(store.attrs["N"])
    task = "INIT"
    for chrom in store.keys():
        counts_dset = store["{}/counts".format(chrom)].value
        missing_rate = counts_dset[:, 3] / float(N)
        missing_rate_dset = store.create_dataset(
            "{}/missing_rates".format(chrom), data=missing_rate)
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
            "HWE": hwe,
            "MISS": missing_rate,
            "AF": af,
            "VAR": var
        }
        # pdb.set_trace()
        msg = pickle.dumps(msg)
        for client in clients:
            Registry.get_instance().set_client_state(client['name'], Commands.INIT_STATS)
            requests.post(f'http://{client["external_host"]}:{client["port"]}/api/init/stats', data=msg)

        hwe_dset = store.create_dataset("{}/hwe".format(chrom), data=hwe)
