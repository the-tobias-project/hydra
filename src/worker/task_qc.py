# stdlib
import pickle
import time

# third party lib
import h5py
import numpy as np


# internal lib
from client.lib import shared
from lib import networking
from lib.utils import write_or_replace
from lib.corr import nancorr
from lib.settings import QCFilterNames, Settings


def init_qc(message, client_config, env):
    print("Pefroming QC")
    filters = pickle.loads(message)
    remove = True
    if "remove" in filters:
        remove = filters["remove"]
    run_QC(filters, client_config, remove=remove, env=env)
    print('Finished reporting counts')


def run_QC(filters, client_config, remove=True, env="production"):
    def find_what_passes(qc_name, dset_name, tokeep, doubleSided=False):
        vals = group[dset_name].value
        if qc_name in filters:
            thresh = float(filters[qc_name])
            if not doubleSided:
                tokeep = np.logical_and(tokeep, vals > thresh)
            else:
                tokeep = np.logical_and(tokeep,
                    np.logical_and(vals > thresh - Settings.kSmallEpsilon,
                        (1.0-vals) > thresh - Settings.kSmallEpsilon))
        return tokeep

    def replace_dataset(tokeep, dset_name, return_deleted=False):
        vals = group[dset_name].value
        remaining = vals[tokeep]
        deleted = vals[np.logical_not(tokeep)]
        write_or_replace(group, dset_name, remaining)
        if return_deleted:
            return deleted
    pfile = client_config["plinkfile"]
    store_name = shared.get_plink_store(pfile)
    with h5py.File(store_name, 'a') as store:
        for chrom in store.keys():
            if chrom == "meta":
                continue
            group = store[chrom]
            positions = group['positions'].value
            tokeep = np.ones_like(positions, dtype=bool)
            tokeep = find_what_passes(QCFilterNames.QC_HWE, "hwe", tokeep)
            tokeep = find_what_passes(QCFilterNames.QC_MAF, "MAF",
                tokeep, doubleSided=True)
            if QCFilterNames.QC_MPS in filters:
                filters[QCFilterName.QC_MPS] = 1 - filters[QCFilterName.QC_MPS]
            tokeep  = find_what_passes(QCFilterNames.QC_MPS, "not_missing_per_snp", tokeep)
            print(f"After filtering {chrom}, {np.sum(tokeep)} snps remain")
            if remove: # Delete what doesn't pass
                replace_dataset(tokeep, 'hwe')
                replace_dataset(tokeep, 'VAR')
                replace_dataset(tokeep, 'MAF')
                replace_dataset(tokeep, 'not_missing_per_snp')
                deleted = replace_dataset(tokeep, 'positions',
                    return_deleted=True)
                for snp in deleted:
                    snp = str(snp)
                    if snp in group:
                        del group[snp]
            else: # Store what has been tagged
                if "passed" in group:
                    tags = group["passed"]
                    tokeep = np.logical_and(tokeep, tags.value)
                    tags[:] = tokeep
                else:
                    if "PCA_passed" in group:
                        del group["PCA_passed"]
                    if "PCA_positions" in group:
                        del group["PCA_positions"]
                    group.create_dataset("PCA_passed",
                        data=np.ones(np.sum(tokeep), dtype=bool))
                    positions = group['positions'].value[tokeep]
                    group.create_dataset("PCA_mask", data=tokeep, dtype=bool)
                    group.create_dataset("PCA_positions", data=positions)
    client_name = client_config['name']
    if remove:
        networking.respond_to_server('api/tasks/QC/FIN', "POST", b'', client_name, env)
    else:
        networking.respond_to_server('api/tasks/PCA/FIN', "POST", b'', client_name, env)


