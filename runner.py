import sys
from subprocess import Popen
import shlex

# Third party lib
from settings import Settings

def main(plinkList, args):
  aux = ""
  for arg in args:
    aux += " {}".format(arg)
  processes = []
  for plinkFile in plinkList:
    arguments = "{} {}".format(plinkFile, aux)
    processes.append(Popen(shlex.split(Settings.python + " client.py " + arguments)))
  for p in processes:
    p.communicate()


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
  main(['testData/dset1', 'testData/dset2', 'testData/dset3'], args)

