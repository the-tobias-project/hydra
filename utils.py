#!/usr/bin/env python3

# stdlib
import _pickle as pickle

# Third party lib
from plinkio import plinkfile
import numpy as np
import h5py
import pdb

# In house lib
from settings import Settings


def encode(message):
    return pickle.dumps(message)


def decode(message):
    return pickle.loads(message)


def write_or_replace(group, name, val, dtype=None):
    if name in group:
        del group[name]
    if not isinstance(val, np.ndarray):
        val = np.array(val)
    if dtype is None:
        dtype = val.dtype
    else:
        val = val.astype(dtype)
    dset = group.create_dataset(name, dtype=dtype, data=val)


def add_pheno(plink_in, multigenecity, out, h=0.85, p_cases=0.5):
    plink_file = plinkfile.open(plink_in)
    if not plink_file.one_locus_per_row():
        print("This script requires that snps are rows and samples columns.")
        exit(1)
    sample_list = plink_file.get_samples()
    locus_list = plink_file.get_loci()
    n = len(sample_list)
    p = len(locus_list)
    edge_offset = 100
    causal_mut_index = np.linspace(
        edge_offset, p-edge_offset, multigenecity, dtype=int)
    gen_effect_size_unnormalized = {item:
        np.random.normal(loc=0, scale=float(h)/np.sqrt(multigenecity))
        for item in causal_mut_index}
    causal_mutations = set()
    mutation_meta = {}
    prs = np.zeros(n)
    for i, variant in enumerate(locus_list):
        row = plink_file.next()
        if i in causal_mut_index:
            genotypes = np.fromiter(row, dtype=float)
            genotypes[genotypes == 3] = np.mean(genotypes[genotypes != 3])
            prs += genotypes * gen_effect_size_unnormalized[i]
    plink_file.close()
    del causal_mut_index, gen_effect_size_unnormalized
    # Draw random environmental effects
    env_rs_unnormalized = np.random.normal(
        loc=0, scale=np.sqrt(1-h**2), size=n)
    gen_effect_size = h * (prs - np.mean(prs)) / np.std(prs)
    env_effect_size = np.sqrt(1-h**2) * (env_rs_unnormalized
        - np.mean(env_rs_unnormalized)) / np.std(env_rs_unnormalized)
    burden = gen_effect_size + env_effect_size
    sorted_i = np.argsort(burden)[::-1]
    ncases = int(n * p_cases)
    cases_i = sorted_i[:ncases]
    affection = np.zeros(n, dtype=np.int8)
    affection[cases_i] = 1
    np.savetxt(out, affection, delimiter='\n', fmt='%1.0i')


def snps_match(plinkName, store_name, position_dset=None):
    # WARNING: this only works if positions are unique.
    plink_file = plinkfile.open(plinkName)
    with h5py.File(store_name, 'r', libver='latest') as store:
        locus_list = plink_file.get_loci()
        plink_file.close()
        plinkSet = set((l.chromosome, l.bp_position) for l in locus_list)
        len_plink = len(plinkSet)
        del locus_list
        if position_dset is None:
            position_dset = 'positions'
        for key in store:
            if key == 'meta':
                continue
            positions = store["{}/{}".format(key, position_dset)].value
            ikey = int(key)
            hset = set((ikey, int(pos)) for pos in positions)
            len_plink -= len(hset)
            plinkSet -= hset
    if len(plinkSet) == 0 and len_plink == 0:
        return True
    return False


def compare_pca(plinkPCA, store_name, dsets_list):
    with h5py.File(store_name, 'r') as store:
        sigmas = store['meta/Sigmas'].value
    with open(plinkPCA+'.eigenval', 'r') as vals_file:
        plinkSigmas = vals_file.read().split()
        plinkSigmas = np.array(plinkSigmas, dtype=float)
    sigmasMatch = np.isclose(sigmas, plinkSigmas, rtol=1e-4, atol=1e-6).all()
    del plinkSigmas, sigmas
    with open(plinkPCA+'.eigenvec', 'r') as vecs_file:
        vecs = vecs_file.read()
        pdb.set_trace()


if __name__ == '__main__':
    np.random.seed(123)
    add_pheno('testData/subsampled', 10, 'testData/wsubsampled.pheno')
