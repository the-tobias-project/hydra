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


def get_covar_file(pfile):
    dirname, basename = os.path.split(pfile)
    prefix = os.path.basename(dirname)
    #print(dirname)
    #current_loc = os.path.dirname(__file__)
    #current_loc = os.path.join("..", "..", current_loc)
    #path = os.path.join(current_loc, dirname)
    #print (f"path is {path}")
    #covar_file = os.path.relpath(current_loc, prefix, "HydraPheno")
    #return covar_file
    path = os.path.join(dirname, "HydraPheno")
    phenopath = os.path.abspath(path)
    return phenopath
