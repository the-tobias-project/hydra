#!/usr/bin/env python3

import argparse
import shlex
from subprocess import Popen
import pdb

# Third party lib
from lib.settings import Settings


names = ["BioME", "MEC_CA", "MEC_HI", "SOL_B", "SOL_M", "SOL_S", "SOL_C", "WHI"]
#names = ["SOL_B", "SOL_M", "SOL_S"]


def worker(plinkList, args):
    aux = ""
    processes = []
    for key in args:
        aux += f'--{key}={args[key]}'
    for i, plinkFile in enumerate(plinkList):
        name = names[i]
        arguments = f"--plinkfile {plinkFile} --name {name} "#{aux}"
        shlex_client = shlex.split(Settings.python + " -m client " + arguments)
        shlex_worker = shlex.split(f"celery -A worker worker -Q {name} -n {name} --concurrency 1")
        processes.append(Popen(shlex_client))
        processes.append(Popen(shlex_worker))
    for p in processes:
        p.communicate()
    return processes


def main():
    parser = argparse.ArgumentParser(description='CWS client')
    parser.add_argument('--local_scratch', type=str, help='Location used for scratch storage during computation')
    parser.add_argument('--port', type=int, help='Which port the clients should use for communication')
    args = parser.parse_args()
    opts = {
        'local_scratch': Settings.local_scratch
    }

    plink_files = [f"../page/splitted/{item}/_ppl" for item in names]
    # plink_files = ['/vagrant/testData/SOL_M/_ppl', '/vagrant/testData/SOL_B/_ppl']
    processes = []
    if args.local_scratch is not None:
        opts['local_scratch'] = args.local_scratch
    if args.port is not None:
        opts['port'] = args.port
    try:
        processes = worker(plink_files, opts)
    except KeyboardInterrupt:
        for p in processes:
            p.kill()
            p.wait(timeout=5)  # Wait up to 5s for process to exit


if __name__ == "__main__":
    main()
