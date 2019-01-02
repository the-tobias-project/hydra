#!/usr/bin/env python3

# std lib
import sys
from subprocess import Popen, PIPE, DEVNULL
import shlex
import time
import os
import shutil
import glob
import tempfile

# Third party lib
from termcolor import colored

# In house Lib
from settings import Settings
from utils import snps_match, compare_pca

import pdb

def process_finished(message):
    if message.startswith("Looks") or message.startswith("Indicate"):
        return True
    return False


def wait_for_process_to_finish(server):
    message = server.stdout.readline()
    while not process_finished(message):
        if message != '':
            print(message)
        message = server.stdout.readline()


def wait_for_client_to_finish(client, k):
    client.stdout.flush()
    for line in client.stdout.readline():
        print(line)
        k -= 1
        client.stdout.flush()


def startup_server_client(scratch=None):
    if scratch is None:
        scratch = Settings.local_scratch
    server = Popen(shlex.split(Settings.python + " server.py " + scratch), stdin=PIPE, 
        stdout=PIPE, bufsize=1, universal_newlines=True)
    client = Popen(shlex.split(Settings.python + " runner.py " + scratch),
        bufsize=1, stdout=PIPE, stderr=DEVNULL, universal_newlines=True)
    message = server.stdout.readline()
    wait_for_process_to_finish(server)
    return server, client


def copy_datasets(location):
    files_to_copy = glob.iglob(os.path.join(Settings.local_scratch, '*.h5py'))
    for f in files_to_copy:
        shutil.copy(f, location)


def test_init():
    server, client = startup_server_client()
    server.stdin.write('init\n')
    print("initalized!")
    wait_for_process_to_finish(server)
    time.sleep(1)
    server.stdin.write('exit\n')
    server.stdin.close()
    return snps_match('testData/subsampled', Settings.local_scratch+'/dset1.h5py')


def run_plink(plink_cmd, inPlink, temp_fldr):
    outname = os.path.join(temp_fldr, os.path.basename(inPlink))
    full_cmd = "{plink} --bfile {plink_file} {cmd} --out {outname}".format(
        plink=Settings.plink, plink_file=inPlink, cmd=plink_cmd, outname=outname)
    print("Running plink command")
    plink_running = Popen(shlex.split(full_cmd), stdin=PIPE, stdout=PIPE)
    plink_running.wait()


def qc_setup(cmd='QC', local_scratch=None):
    if local_scratch is None:
        local_scratch = Settings.local_scratch
    temp_location = tempfile.mkdtemp(prefix=local_scratch+"/")
    copy_datasets(temp_location)
    server, client = startup_server_client(temp_location)
    server.stdin.write('{}\n'.format(cmd))
    return temp_location, server, client


def test_qc_hwe(threshold):
    temp_location, server, client = qc_setup()
    server.stdin.write('hwe {}\n'.format(threshold))
    wait_for_process_to_finish(server)
    plink_cmd = "--hwe {} midp --make-bed".format(threshold)
    run_plink(plink_cmd, 'testData/subsampled', temp_location)
    server.stdin.write('exit\n')
    server.stdin.close()
    plink_to_compare_to = os.path.join(temp_location, 'subsampled')
    results = snps_match(plink_to_compare_to, temp_location+'/central.h5py')
    shutil.rmtree(temp_location)
    return results


def test_qc_maf(threshold):
    temp_location, server, client = qc_setup()
    server.stdin.write('maf {}\n'.format(threshold))
    wait_for_process_to_finish(server)
    plink_cmd = "--maf {} --make-bed".format(threshold)
    run_plink(plink_cmd, 'testData/subsampled', temp_location)
    server.stdin.write('exit\n')
    server.stdin.close()
    time.sleep(.1)
    plink_to_compare_to = os.path.join(temp_location, 'subsampled')
    results = snps_match(plink_to_compare_to, temp_location+'/central.h5py')
    shutil.rmtree(temp_location)
    return results


def test_qc_mps(threshold):
    temp_location, server, client = qc_setup()
    server.stdin.write('mps {}\n'.format(threshold))
    wait_for_process_to_finish(server)
    plink_cmd = "--geno {} --make-bed".format(threshold)
    run_plink(plink_cmd, 'testData/subsampled', temp_location)
    server.stdin.write('exit\n')
    server.stdin.close()
    time.sleep(.1)
    plink_to_compare_to = os.path.join(temp_location, 'subsampled')
    results = snps_match(plink_to_compare_to, temp_location+'/central.h5py')
    shutil.rmtree(temp_location)
    return results


##### PCA TESTS


def test_pca_ld_pruning(win, num_pcs):
    temp_location, server, client = qc_setup('pca')
    server.stdin.write('maf 0.1 hwe 1e-5 ld {}\n'.format(win))
    wait_for_process_to_finish(server)
    server.stdin.write('exit\n')
    server.stdin.close()
    plink_cmd = "--maf 0.1 --hwe 1e-5 midp --indep-pairwise {} 25 0.2".format(win)
    run_plink(plink_cmd, 'testData/subsampled', temp_location)
    plink_cmd = "--extract {}/subsampled.prune.in --make-bed".format(
        temp_location)
    run_plink(plink_cmd, 'testData/subsampled', temp_location)
    plink_to_compare_to = os.path.join(temp_location, 'subsampled')
    ld_results = snps_match(plink_to_compare_to, temp_location+'/central.h5py', 
        'PCA_positions')
    ## Now we check theactual pcs
    plink_cmd = "--pca {}".format(num_pcs)
    plink_loc = temp_location+'/subsampled'
    run_plink(plink_cmd, temp_location+'/subsampled', temp_location)
    dsets = [temp_location+'/dset1.h5py', temp_location+'/dset2.hpy', 
        temp_location+'/dset3.h5py']
    compare_pca(plink_loc, temp_location+'/central.h5py', dsets)
    return ld_results, temp_location



def run_tests():
    assert test_init(), "Initialization failed"
    print(colored("Initialization test: ",'red'), colored(u'\u2713', 'red'))
    assert test_qc_hwe(1e-5), "HWE failed"
    print(colored("QC HWE test: ",'red'), colored(u'\u2713', 'red'))
    assert test_qc_maf(0.05), "MAF failed"
    print(colored("QC maf test: ",'red'), colored(u'\u2713', 'red'))
    assert test_qc_mps(0.05), "Missing per snp failed"
    print(colored("QC missing per snp test: ",'red'), colored(u'\u2713', 'red'))
    results, pca_temp_location = test_pca_ld_pruning(50, 10)
    assert results, "LD pruning failed"
    print(colored("LD pruning test: ",'red'), colored(u'\u2713', 'red'))


def main():
    run_tests()


if __name__ == '__main__':
    main()

