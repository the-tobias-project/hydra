#!/usr/bin/env python3

import argparse
import random
import shlex
<<<<<<< HEAD
from subprocess import Popen
=======
import pdb
>>>>>>> pre-switch to http

# Third party lib
from src.lib.settings import Settings


names = ['zealot', 'dragoon', 'scout', 'disruptor', 'probe', 'nexus',
         'pylon', 'carrier', 'arbiter', 'immortal', 'colossus', 'templar', 'archon']
processes = []

<<<<<<< HEAD
=======
if __name__ == "__main__":
  if len(sys.argv) >= 2:
    local_scratch = sys.argv[1]
    if len(sys.argv) == 3:
        PORT = int(sys.argv[2])
  else:
    local_scratch = Settings.local_scratch
  args = [local_scratch]
  if len(sys.argv) == 3:
    args += [PORT] 
  #main(['testData/dset1', 'testData/dset2', 'testData/dset3'], args)
  main(['page/splitted/BioME/_ppl', "page/splitted/MEC_CA/_ppl",
    'page/splitted/MEC_HI/_ppl', 'page/splitted/SOL_B/_ppl'
    , 'page/splitted/SOL_C/_ppl', 'page/splitted/SOL_M/_ppl', 
    'page/splitted/SOL_S/_ppl', 'page/splitted/WHI/_ppl'], args)
>>>>>>> pre-switch to http

def worker(plinkList, args):
    aux = ""
    for key in args:
        aux += f'--{key}={args[key]}'
    for plinkFile in plinkList:
        name_index = random.randint(0, len(names) - 1)
        name = names.pop(name_index)
        arguments = f"{plinkFile} {aux} --name={name}"
        shlex_split = shlex.split(Settings.python + " client.py " + arguments)
        print(shlex_split)
        global processes
        processes.append(Popen(shlex_split))
    for p in processes:
        p.communicate()


def main():
    parser = argparse.ArgumentParser(description='CWS client')
    parser.add_argument('--local_scratch', type=str, help='Location used for scratch storage during computation')
    parser.add_argument('--port', type=int, help='Which port the clients should use for communication')
    args = parser.parse_args()
    opts = {
        'local_scratch': Settings.local_scratch
    }

    # files = ['testData/popres1', 'testData/popres2', 'testData/popres3']
    files = ['testData/popres1']
    if args.local_scratch is not None:
        opts['local_scratch'] = args.local_scratch
    if args.port is not None:
        opts['port'] = args.port
    try:
        worker(files, opts)
    except KeyboardInterrupt:
        for p in processes:
            p.kill()


if __name__ == "__main__":
    main()
