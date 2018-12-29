import sys
from subprocess import Popen
import shlex

# Third party lib
from settings import Settings

def main(plinkList, local_scratch):
  processes = []
  for plinkFile in plinkList:
    arguments = "{} {}".format(plinkFile, local_scratch)
    processes.append(Popen(shlex.split(Settings.python + " client.py " + arguments)))
  for p in processes:
    p.communicate()


if __name__ == "__main__":
  if len(sys.argv) == 2:
    local_scratch = sys.argv[1]
  else:
    local_scratch = Settings.local_scratch
  main(['testData/dset1', 'testData/dset2', 'testData/dset3'], local_scratch)

