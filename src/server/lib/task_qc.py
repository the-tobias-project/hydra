# stdlib
import logging
import pickle
import re
import os
import pdb

# third party lib
import requests
import h5py
import numpy as np

# internal lib
from lib.settings import Settings, Options, Commands
from lib.client_registry import Registry


storePath = os.path.join(Settings.local_scratch, "central.h5py")
store = h5py.File(storePath, "a")
clients = Registry.get_instance().list_clients()


def split_command(command):
    command = command.upper()
    refloat = "[0-9,E,.,-]*"
    filters = {}
    name = Options.HWE
    x = re.search(name+refloat, command)
    if x: # hwe filter
        filters[name] = float(x.group()[len(name):])
    name = Options.MAF
    x = re.search(name+refloat, command)
    if x: # maf filter
        filters[name] = float(x.group()[len(name):])
    name = Options.MPS
    x = re.search(name+refloat, command)
    if x: # missing per SNP is filter
        filters[name] = float(x.group()[len(name):])
    name = Options.LD
    x = re.search(name+"[0-9]*_"+refloat, command)
    if x:
        x = x.group()[len(name):]
        filters[name] = x.split("_")
    return filters


def start_client_qc_task(filters, stage=Commands.QC):
    for client in clients:
        Registry.get_instance().set_client_state(client['name'], stage)
        data = pickle.dumps(filters)
        requests.post(f'http://{client["external_host"]}:{client["port"]}/api/qc', data=data)


def start_local_qc_task(filters, prefix=None): #Filter based on local info
    """Performs the filters with threshold values specified in a dictionary
    named filters. This deletes the snps if prefix is left as None"""
    for chrom in store.keys():
        if chrom == 'meta':
            continue
        group  = store[chrom]
        pos    = group['positions']
        counts = group['counts']
        mr     = group['missing_rates']
        af     = group['allele_freq']
        hwe    = group['hwe']
        tokeep = np.ones(shape=pos.value.shape, dtype=bool)
        if Options.HWE in filters:
            val = filters[Options.HWE]
            tokeep = np.logical_and(tokeep, hwe.value > val)
        if Options.MAF in filters:
            val = filters[Options.MAF]
            tokeep = np.logical_and(tokeep, af.value > val - Settings.kSmallEpsilon)
            tokeep = np.logical_and(tokeep, 1.0-af.value > val - Settings.kSmallEpsilon)
        if Options.MPS in filters:
            val = filters[Options.MPS]
            tokeep = np.logical_and(tokeep, mr.value < val)
        logging.info("In chromosome {}, {} snps were deleted and {} snps remain".format(chrom,
            tokeep.shape[0] - np.sum(tokeep), np.sum(tokeep)))
        # Delete or tag the filtered locations
        if prefix is None:
            pos_vals, counts_vals = pos.value[tokeep], counts.value[tokeep]
            mr_vals = mr.value[tokeep]
            del group["positions"], group["counts"], group["missing_rates"]
            d1 = group.require_dataset("positions", pos_vals.shape
                , dtype = pos_vals.dtype)
            d1[:] = pos_vals
            d2 = group.require_dataset("counts", counts_vals.shape
                , dtype=counts_vals.dtype)
            d2[:] = counts_vals
            d3 = group.require_dataset("missing_rates", mr_vals.shape
                , dtype=mr_vals.dtype)
            d3[:] = mr_vals
            del pos_vals, counts_vals, mr_vals
            af_vals, hwe_vals= af.value[tokeep], hwe.value[tokeep]
            del group["hwe"], group["allele_freq"]
            d4 = group.require_dataset("hwe", hwe_vals.shape
                , dtype=hwe_vals.dtype)
            d4[:] = hwe_vals
            d5 = group.require_dataset("allele_freq", af_vals.shape
                , dtype=af_vals.dtype)
            d5[:] = af_vals
        else:
            n = np.sum(tokeep)
            ones = np.ones(n, dtype=bool)
            d1 = group.require_dataset(prefix + "passed", ones.shape
                , dtype=bool)
            d1[:] = ones
            pos_vals = pos.value[tokeep]
            d2 = group.require_dataset(prefix + "positions", pos_vals.shape
                , dtype=pos_vals.dtype)
            d2[:] = pos_vals
            af_vals = af.value[tokeep]
            d3 = group.require_dataset(prefix + "allele_freq"
                , af_vals.shape, dtype=af_vals.dtype)
            d3[:] = af_vals

def filter_finished(client_name, state):
    Registry.get_instance().set_client_state(client_name, "Filterd")
    if not Registry.get_instance().num_clients_in_state(state):
        logging.info(f"Done with filtering in {Commands.QC} stage")
        return True
    return False
