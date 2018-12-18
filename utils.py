#!/usr/bin/env python3

import _pickle as pickle

import numpy as np


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
