# stdlib
import pickle
import logging

# third party lib
import h5py
import numpy as np


# internal lib
from client.lib import shared
from lib import networking
from lib.utils import write_or_replace
from lib.settings import QCFilterNames, Settings

logger = logging.getLogger("worker")


def init_qc(message, client_config, env):
    logger.info("Pefroming QC.")
    filters = pickle.loads(message)
    remove = filters.get("remove", False)
    run_QC(filters, client_config, remove=remove, prefix=filters.get("mask_prefix", None), env=env)
    logger.info('Finished reporting counts.')


def run_QC(filters, client_config, prefix, remove=True, env="production"):
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
                filters[QCFilterNames.QC_MPS] = 1 - filters[QCFilterNames.QC_MPS]
            tokeep = find_what_passes(QCFilterNames.QC_MPS, "not_missing_per_snp", tokeep)
            logger.info(f"After filtering {chrom}, {np.sum(tokeep)} snps remain")
            if remove:  # Delete what doesn't pass
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
            else:  # Store what has been tagged
                pass_mask = prefix + "_mask"
                pos_mask = prefix + "_positions"
                if pass_mask in group:
                    del group[pass_mask]
                if pos_mask in group:
                    del group[pos_mask]
                write_or_replace(group, pass_mask, val=tokeep, dtype=bool)
                positions = group['positions'].value[tokeep]
                write_or_replace(group, pos_mask, val=positions)
                if prefix == "PCA":
                    write_or_replace(group, "PCA_passed", val=np.ones(np.sum(tokeep), dtype=bool))
                    if 'non_ld_mask' in group:
                        del group['non_ld_mask']
    client_name = client_config['name']
    if prefix == "QC":
        networking.respond_to_server('api/tasks/QC/FIN', "POST", b'', client_name, env)
    else:
        networking.respond_to_server('api/tasks/PCA/FIN', "POST", b'', client_name, env)
