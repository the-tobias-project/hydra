import sys
from subprocess import Popen
import shlex
import pdb

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
  #main(['testData/dset1', 'testData/dset2', 'testData/dset3'], args)
  main(['page/splitted/BioME/_ppl', "page/splitted/MEC_CA/_ppl",
    'page/splitted/MEC_HI/_ppl', 'page/splitted/SOL_B/_ppl'
    , 'page/splitted/SOL_C/_ppl', 'page/splitted/SOL_M/_ppl', 
    'page/splitted/SOL_S/_ppl', 'page/splitted/WHI/_ppl'], args)

