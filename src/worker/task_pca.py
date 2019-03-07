# stdlib
import os
import pickle
import socket
import sys
import time

# third party lib
from celery import current_app
import h5py
import numpy as np
from plinkio import plinkfile


# internal lib
from client.lib import shared
from lib import HTTPResponse
from lib.utils import write_or_replace
from lib.corr import nancorr, process_plink_row


def report_ld(data, client_config):
    pfile = client_config["plinkfile"]
    with h5py.File(shared.get_plink_store(pfile), 'a') as store:
        n = store.attrs['n']
        if self.r3 == 0: # fake start the first round
            for chrom in self.chroms: 
                # if the length is less than r1, you deserve an error. 
                # No apologies 
                tags = store["{}/PCA_passed".format(chrom)]
                message[chrom] = tags[0:self.r1]
        for key, state in message.items():
            if key == "TASK" or key == "SUBTASK":
                continue
            chrom = key
            tags = store["{}/PCA_passed".format(chrom)]
            if state[0] == "E": # Finished with this chrom
                if len(message) == 3: # Done with everything
                    self.do_ld = False
                    self.chroms = None
                    msg = {"TASK": "PCA", "SUBTASK": "PCA_POS"}
                    self.server.message(encode(msg))
                    self.verbose = True
                    return
                continue
            else:
                tokeep = state
                end = self.r3 + len(tokeep)
            pos = store["{}/PCA_positions".format(chrom)]
            positions = pos[self.r3: end]
            positions = positions[tokeep]
            genotypes = np.empty((n,len(positions)), dtype=np.float32)
            for i, snp in enumerate(positions):
                genotypes[:,i] = store["{}/{}".format(chrom, snp)].value
            corr = nancorr(genotypes)
            msg[chrom] = corr
    self.server.message(encode(msg))
    self.r3 += self.r2        

def init_store(client_config):
    plinkToH5(client_config)
    print("preparing counts")
    report_counts(client_config)
    print('Finished reporting counts')


def send_counts_to_server(data, client_config):
    print('sending counts to server')
    client_name = client_config['name']
    data = pickle.dumps(data)
    HTTPResponse.respond_to_server('api/tasks/INIT/COUNT', 'POST', data, client_name)
    print('sent counts to server')


def init_stats(message, client_config):
    print('Inside init_stats')
    # Wait on previous tasks to finish
    i = current_app.control.inspect()
    client_name = client_config['name']
    while True:
        active_tasks = i.active()[f'celery@{client_name}']
        dependent_tasks = list(filter(lambda x: x['name'] == 'tasks.init_store', active_tasks))
        if len(dependent_tasks) > 0:
            print('Waiting on tasks.init_store to finish')
            time.sleep(.1)
        else:
            break
    print('Resuming with init_stats')
    message = pickle.loads(message)
    chrom = message["CHROM"]
    print(f'Chrom: {chrom}')
    pfile = client_config['plinkfile']
    with h5py.File(shared.get_plink_store(pfile), 'a') as store:
        chrom_group = store[str(chrom)]
        if "MISS" in message:
            vals = message["MISS"]
            task = "not_missing_per_snp"
            dset = chrom_group.create_dataset(task, data=1 - vals)
        if "AF" in message:
            vals = message["AF"]
            task = 'MAF'
            dset = chrom_group.create_dataset(task, data=vals)
        if "HWE" in message:
            vals = message["HWE"]
            task = "hwe"
            dset = chrom_group.create_dataset(task, data=vals)
        if "VAR" in message:
            vals = message["VAR"]
            task = "VAR"
            dset = chrom_group.create_dataset(task, data=vals)
    print('Finished with init_stats')

    client_name = client_config['name']
    status = f'Finished with init stats for chrom {chrom}'
    HTTPResponse.respond_to_server(f'api/clients/{client_name}/report?status={status}', 'POST')
