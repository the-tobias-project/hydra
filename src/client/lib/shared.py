# stdlib
import logging
import os

# internal lib
from lib.settings import Settings


# -=-=-=-=-= Plink-related functions -=-=-=-=-=-
def set_plinkfile(pfile):
    global plinkfile
    plinkfile = pfile


def get_store_path(plink_name):
    dirname, basename = os.path.split(plink_name)
    prefix = os.path.basename(dirname)
    write_to = Settings.local_scratch
    store_name = os.path.join(write_to, prefix + basename + '.h5py')
    return store_name


def get_plink_store(pfile):
    dirname, basename = os.path.split(pfile)
    prefix = os.path.basename(dirname)
    write_to = Settings.local_scratch
    store_name = os.path.join(write_to, prefix + basename + '.h5py')
    return store_name


def get_covar_file():
    dirname, basename = os.path.split(plinkfile)
    prefix = os.path.basename(dirname)
    covar_file = os.path.join(prefix, "HydraPheno")
    return covar_file
