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
    return pickle.dumps(message) + Settings.PICKLE_STOP


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
    fids = np.array([item.fid for item in sample_list])
    iids = np.array([item.iid for item in sample_list])
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
    affection[cases_i] = 2
    affection[affection == 0] = 1
    towrite = np.column_stack((fids, iids, affection))
    np.savetxt(out, towrite, delimiter='\t', fmt=['%s', '%s', '%s'], header='FID\tID\tpheno',)


def snps_match(plinkName, store_name, position_dset=None):
    # WARNING: this only works if positions are unique.
    with h5py.File(store_name, 'r', libver='latest') as store:
        # check the plink file
        plink_file = plinkfile.open(plinkName)
        locus_list = plink_file.get_loci()
        plink_file.close()
        plinkSet = set((l.chromosome, l.bp_position) for l in locus_list)
        del locus_list
        len_plink = len(plinkSet)
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
        k = len(sigmas)
    sigmasMatch = np.isclose(sigmas, plinkSigmas, rtol=1e-4, atol=1e-6).all()
    del plinkSigmas, sigmas
    plinkVals = np.loadtxt(plinkPCA+'.eigenvec', usecols=range(2,2+k))
    valsMatch = True
    index_new, index_old = 0, 0
    for dset in dsets_list:
        with h5py.File(dset, 'r') as store:
            vals = store['meta/pca_u'].value
            index_new += vals.shape[0]
            valsMatch * np.isclose(plinkVals[index_old:index_new,:], vals, rtol=1e-4, atol=1e-6).all()
            index_old = index_new
    return valsMatch * sigmasMatch


def compare_regression(plinkRegression, store_name, dset_name="results"):
    converter = lambda x: np.nan if x == b'NA' else float(x)
    plinkResults = np.loadtxt(plinkRegression, usecols=[0, 6,7,8], skiprows=1, 
        converters={6: converter, 7: converter, 8: converter})
    with h5py.File(store_name, 'r') as store:
        pdb.set_trace()
        for key in store:
            if key != "meta":
                results = store[key+"/{}".format(dset_name)].value
                res = np.array([i for subset in results for i in subset])
                plinkSubset_ind = plinkResults[:,0] == int(key)
                plinkSubset = plinkResults[plinkSubset_ind, :]
                del plinkSubset_ind
                
          

    pass
    


if __name__ == '__main__':
    np.random.seed(123)
    add_pheno('testData/subsampled', 10, 'testData/subsampled.pheno')
