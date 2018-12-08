import sys
from subprocess import Popen
import shlex
PYTHON="/srv/gsfs0/software/python/3.6.4/bin/python3"

def main(plinkList):
  processes = []
  for plinkFile in plinkList:
    arguments = plinkFile 
    print (arguments)
    processes.append(Popen(shlex.split(PYTHON + " client.py " + arguments)))
  for p in processes:
    p.wait()


if __name__ == "__main__":
  main(['smallData/popres1', 'smallData/popres2', 'smallData/popres3'])

