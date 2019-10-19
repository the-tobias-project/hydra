# stdlib
import os
import pickle
import sys
import time
import subprocess
import logging

# third party lib
from celery import current_app
import h5py
import numpy as np
from plinkio import plinkfile

# internal lib
from client.lib import shared
from lib import networking
from lib.utils import write_or_replace
from lib.corr import process_plink_row

logger = logging.getLogger("worker")


def init_store(client_config, env):
    pfile = client_config['plinkfile']
    store_name = shared.get_plink_store(pfile)
    if os.path.isfile(store_name):
        clear_consistency_flag(store_name)
        report_file_info(store_name, client_config, env)
        logger.info(f"HDF5 file {store_name} already exists.")
    else:
        plinkToH5(client_config, env)
    logger.info("Preparing counts.")
    report_counts(client_config, env)
    logger.info('Finished reporting counts.')


def clear_consistency_flag(fname):
    try:
        child_sig = subprocess.call([f"h5clear -s {fname}"], shell=True)
        if child_sig != 0:
            child_sig("Consistency flag has not been cleared.", child_sig, file=sys.stderr)
    except OSError as e:
        logger.Error("Execution failed:", e, file=sys.stderr)


def report_file_info(store_name, client_config, env):
    # TODO make QC not actually delete stuff
    with h5py.File(store_name, 'a') as store:
        store.attrs['has_global_AF'] = False
        store.attrs['has_centering'] = False
        store.attrs['has_normalization'] = False
        if "meta" in store:
            del store["meta"]
        if "pca" in store:
            del store["pca"]
        for chrom in store:
            positions = store[f"{chrom}/positions"].value
            send_positions_to_server(positions, chrom, client_config, env)


def plinkToH5(client_config, env):
    """Gets plink prefix, produces an HDF file with the same prefix"""
    pfile = client_config['plinkfile']
    store_name = shared.get_plink_store(pfile)
    logger.info(f'Opening plinkfile: {pfile}')
    try:
        plink_file = plinkfile.open(pfile)
    except MemoryError as e:
        logger.error('MemoryError!')
        logger.error(e)
    if not plink_file.one_locus_per_row():
        logger.error("""This script requires that snps are
            rows and samples columns.""")
        sys.exit(1)
    sample_list = plink_file.get_samples()
    locus_list = plink_file.get_loci()
    n_tot = len(sample_list)
    logger.info(f'Opening h5py file:{store_name}')
    with h5py.File(store_name, 'w', libver='latest') as store:
        store.attrs['n'] = len(sample_list)
        store.attrs['has_local_AF'] = False
        store.attrs['has_global_AF'] = False
        store.attrs['has_centering'] = False
        store.attrs['has_normalization'] = False
        potential_pheno_file = pfile+".pheno"
        if os.path.isfile(pfile+".pheno"):
            affection = np.loadtxt(potential_pheno_file, dtype=int, usecols=2)
        else:
            affection = [sample.affection for sample in sample_list]
        if len(np.unique(affection)) > 2:
            raise ValueError("phenotype is not binary. We only support binary for now")
        write_or_replace(store, 'meta/Status', affection, np.int8)
        ids = [sample.iid for sample in sample_list]
        write_or_replace(store, 'meta/id', ids, 'S11')
        del ids, affection
        # Read Demographic file
        logger.info(f'Reading demographic file at {pfile}.ind')
        logger.info(f'File exists: {os.path.isfile(pfile + ".ind")}')
        with open(pfile + ".ind", 'r') as dem_f:
            dem = [(row.split("\t")[2]).encode("UTF8") for row in dem_f]
            write_or_replace(store, 'meta/regions', dem)
        # Read chromosome data
        current_chr = 1
        positions = []
        rsids = []
        all_counts = []
        current_group = store.require_group(str(current_chr))
        genotypes = np.zeros(n_tot, dtype=np.float32)
        for locus, row in zip(locus_list, plink_file):
            if locus.chromosome != current_chr:
                if len(positions) == 0:
                    del store[str(current_chr)]
                else:
                    write_or_replace(current_group, 'positions', positions,
                                     dtype=np.uint)
                    write_or_replace(current_group, 'rsids', rsids)
                    write_or_replace(current_group, 'counts', all_counts,
                                     np.uint32)

                    send_positions_to_server(positions, current_chr, client_config, env)
                    positions = []
                    # rsid = []
                    all_counts = []
                current_chr = locus.chromosome
                if current_chr == 23:
                    break
                current_group = store.require_group(str(current_chr))
            pos = str(locus.bp_position)
            counts, geno = process_plink_row(row, genotypes)
            # This should be a try except
            try:
                current_group.create_dataset(pos, data=geno)
            except Exception:
                logger.error(f"Cannot write position: chr{locus.chromosome} {pos}")
            rsids.append(locus.name.encode('utf8'))
            positions.append(pos)
            all_counts.append(counts)
        if locus.chromosome != 23:
            write_or_replace(current_group, 'positions', positions,
                             np.uint32)
            write_or_replace(current_group, 'rsids', rsids)
            write_or_replace(current_group, 'counts', all_counts, np.uint32)
            send_positions_to_server(positions, current_chr, client_config, env)
    plink_file.close()
    logger.info('Finished writing plink to hdf5.')


def report_counts(client_config, env):
    """
    Report the counts (Het, homo Alt, missing)
    """
    pfile = client_config['plinkfile']
    store_name = shared.get_plink_store(pfile)
    with h5py.File(store_name, 'r') as store:
        countDict = {}
        n = store.attrs["n"]
        countDict["START"] = True
        keys = [i for i in store.keys() if i != 'meta']
        for chrom in keys:
            countDict["n"] = int(n)
            countDict["CHROM"] = chrom
            count_arr = store["{}/counts".format(chrom)].value
            countDict["COUNTS"] = count_arr
            if chrom == keys[-1]:
                countDict["END"] = True
            logger.info(f'Sending counts from chrom {chrom}.')
            send_counts_to_server(countDict, client_config, env)


def send_positions_to_server(positions, chrom, client_config, env):
    client_name = client_config['name']

    data = pickle.dumps({
        'CHROM': chrom,
        'POS': positions
    })

    networking.respond_to_server('api/tasks/INIT/POS', 'POST', data, client_name, env)


def send_counts_to_server(data, client_config, env):
    client_name = client_config['name']
    data = pickle.dumps(data)
    networking.respond_to_server('api/tasks/INIT/COUNT', 'POST', data, client_name, env)


def init_stats(message, client_config, env):
    # Wait on previous tasks to finish
    i = current_app.control.inspect()
    client_name = client_config['name']
    while i.active() is not None:
        active_tasks = i.active()[f'celery@{client_name}']
        dependent_tasks = list(filter(lambda x: x['name'] == 'tasks.init_store', active_tasks))
        if len(dependent_tasks) > 0:
            logger.info('Waiting on tasks.init_store to finish.')
            time.sleep(.1)
        else:
            break
    message = pickle.loads(message)
    chrom = message["CHROM"]
    logger.info(f'Computing statistics for Chrom: {chrom}.')
    pfile = client_config['plinkfile']
    with h5py.File(shared.get_plink_store(pfile), 'a') as store:
        chrom_group = store[str(chrom)]
        if "MISS" in message:
            vals = message["MISS"]
            task = "not_missing_per_snp"
            write_or_replace(chrom_group, task, val=1-vals)
        if "AF" in message:
            vals = message["AF"]
            task = 'MAF'
            write_or_replace(chrom_group, task, val=vals)
        if "HWE" in message:
            vals = message["HWE"]
            task = "hwe"
            write_or_replace(chrom_group, task, val=vals)
        if "VAR" in message:
            vals = message["VAR"]
            task = "VAR"
            write_or_replace(chrom_group, task, val=vals)
    logging.info(f'Finished initializing QC statistics for chrom {chrom}.')

    client_name = client_config['name']
    status = f'Finished with init stats for chrom {chrom}.'
    networking.respond_to_server(f'api/clients/{client_name}/report?status={status}', 'POST', env=env)
