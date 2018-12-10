# Armin Pourshafeie

#TODO documentation
# Running the gwas
# Careful here, eigh uses https://software.intel.com/en-us/mkl-developer-reference-c-syevr behind the hood
# so it can be significantly slower 



import _pickle as pickle
import pdb
import tqdm, time, json
import numpy as np
import h5py, os, sys
from plinkio import plinkfile

# In house packages files
from corr import nancorr, corr, hweP
#from optimizationAux import * # This will be used later

with open('GLOBALS.json', 'r') as fp:
  locals_dict = json.load(fp)
for key,val in locals_dict.items():
    exec(key + '=val')


kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;
kSmallEpsilon = 0.00000000000005684341886080801486968994140625;
kLargeEpsilon = 1e-7
SMALL_EPSILON = 0.00000000000005684341886080801486968994140625


def encode(data):
  return pickle.dumps(data)

def decode(data):
  return pickle.loads(data)

class ServerTalker(object):
  """Dispatches commands from/to the server"""
  def __init__(self, connections, server, scratch):
    self.connections = connections
    self.status      = "Initialized"
    self.task        = None
    self.subtasks    = None
    self.server      = server
    self.scratch     = scratch
    self.storePath   = os.path.join(scratch, "central.h5py")
    self.store       = h5py.File(self.storePath, "a")
    self.counter     = connections
    self.r0          = None
    self.r1, self.r2 = None, None
    self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
    self.tqdm        = None

  def report_status(self):
    print("Now working on: {}".format(self.status))
  def dispatch(self,message):
    task = message["TASK"]
    self.task = task
    subtask = message["SUBTASK"]
    self.subtask = subtask
    current_task = "{} {}".format(task, subtask)
    if self.status != current_task:
      self.status = current_task
      self.report_status()
    if self.task == "INIT":
      if subtask == "POS":
        self.store_positions(message)
      elif subtask == "COUNT":
        self.store_counts(message)
        print("Count statistics have been initialized!")
    elif self.task == "QC":
      self.QC_filters(message["SUBTASK"], message["VALS"])
    elif self.task == "PCA":
      if self.subtask == "FILTERS":
        self.PCA_filters(message["FILTERS"], message["VALS"])
        self.counter = self.connections
        self.r0      = 0
        return
      if self.subtask == "LD":
        self.ld_filters(message)
      if self.subtask == "PCA_POS":
        self.counter -= 1
        if self.counter == 0:
          msg = {"TASK":self.task, "SUBTASK":self.subtask}
          self.report_content("PCA_passed", msg)
      if self.subtask == "COV":
        self.buildCov(message)#START HERE

    elif task == "Association":
      pass 

  def store_positions(self, message):
    chrom = message["CHROM"] 
    positions = message["POS"]
    dsetname = "{}/positions".format(chrom)
    if dsetname not in self.store:
      print("{} loci in chromosome {}".format(len(positions), chrom))
      dset = self.store.require_dataset(dsetname, (len(positions),), 
          dtype=np.uint32)
      dset[:] = np.array(positions, dtype=np.uint32)


  def store_counts(self, message):
    n = message["n"] 
    if "START" in message:
      if "N" not in self.store.attrs:
        self.store.attrs["N"] = 0
      self.store.attrs["N"] += n
    chrom = message["CHROM"]
    size = len(message["COUNTS"])
    dsetname = "{}/counts".format(chrom)
  
    if dsetname not in self.store:
      dset = self.store.require_dataset(dsetname, (size, 4), dtype=np.uint32)
    else :
      dset = self.store[dsetname]
    counts = message["COUNTS"]
    homo_ref = n - np.sum(counts, axis = 1)[:, np.newaxis].astype(np.uint32)
    dset[:] += np.hstack((homo_ref, counts))
    if "END" in message:
      self.counter -= 1
      if self.counter == 0: # Everybody's report is in
        self.status = "INIT STATS"
        self.report_status()
        self.count_stats()
        self.server.run_QC()
  
  def count_stats(self):
    N = float(self.store.attrs["N"])
    task = "INIT"

    for chrom in self.store.keys():
      counts_dset = self.store["{}/counts".format(chrom)].value
      missing_rate = counts_dset[:,3] / float(N)
      #msg = {"TASK":task, "SUBTASK": "STATS", "CHROM": chrom, "MISS": missing_rate}
      #self.server.message(encode(msg))
      missing_rate_dset = self.store.create_dataset("{}/missing_rates".format(chrom), data=missing_rate)
      af = (counts_dset[:,2] * 2 + counts_dset[:,1]).astype(float)
      af /= (np.sum(counts_dset[:,:3], axis=1)*2)
      af = af #np.minimum(af, 1-af)
      #msg = {"TASK":task, "SUBTASK": "STATS", "CHROM":chrom, "AF": af}
      #self.server.message(encode(msg))
      self.store.create_dataset("{}/allele_freq".format(chrom), data=af)
      var = counts_dset[:,0] * (2*af)**2 + counts_dset[:,1] * (1-2*af)**2 + counts_dset[:,2] * (2-2*af)**2
      var /= (N-counts_dset[:,3]) # 2*af*(1-af)
      self.store.create_dataset("{}/var".format(chrom), data=var)
      hwe = hweP(counts_dset[:,:3].astype(np.int32), 1) # Recompile HWEP with uint32
      msg = {"TASK":task, "SUBTASK": "STATS", "CHROM":chrom, "HWE": hwe, "MISS":missing_rate, "AF":af, "VAR":var}
      self.server.message(encode(msg))
      hwe_dset = self.store.create_dataset("{}/hwe".format(chrom), data=hwe)



### QC


  def QC_filters(self, filter_list, value_list, prefix=None):
    for chrom in self.store.keys():
      group  = self.store[chrom]
      pos    = group['positions']
      counts = group['counts']
      mr     = group['missing_rates']
      af     = group['allele_freq']
      hwe    = group['hwe']
      tokeep = np.ones(shape=pos.value.shape, dtype=bool)
      if QC_HWE in filter_list:
        ind = filter_list.index(QC_HWE)
        tokeep = np.logical_and(tokeep, hwe.value > value_list[ind])
      if QC_MAF in filter_list:
        ind = filter_list.index(QC_MAF)
        tokeep = np.logical_and(tokeep, af.value > value_list[ind])
        tokeep = np.logical_and(tokeep, 1.0-af.value > value_list[ind])
      if QC_MPS in filter_list:
        ind = filter_list.index(QC_MPS)
        tokeep = np.logical_and(tokeep, mr.value < value_list[ind])
      print("in chromosome {}, {} snps were deleted and {} snps remain".format(chrom, tokeep.shape[0] - np.sum(tokeep), np.sum(tokeep)))
        
      if prefix is None: 
        pos_vals, counts_vals, mr_vals = pos.value[tokeep], counts.value[tokeep], mr.value[tokeep]
        del group["positions"], group["counts"], group["missing_rates"]
        d1 = group.require_dataset("positions", pos_vals.shape, dtype = pos_vals.dtype)
        d1[:] = pos_vals
        d2 = group.require_dataset("counts", counts_vals.shape, dtype=counts_vals.dtype)
        d2[:] = counts_vals
        d3 = group.require_dataset("missing_rates", mr_vals.shape, dtype=mr_vals.dtype)
        d3[:] = mr_vals
        del pos_vals, counts_vals, mr_vals
        af_vals, hwe_vals= af.value[tokeep], hwe.value[tokeep]
        del group["hwe"], group["allele_freq"]
        d4 = group.require_dataset("hwe", hwe_vals.shape, dtype=hwe_vals.dtype)
        d4[:] = hwe_vals
        d5 = group.require_dataset("allele_freq",af_vals.shape, dtype=af_vals.dtype)
        d5[:] = af_vals

      else:
        d1 = group.require_dataset(prefix + "passed", tokeep.shape, dtype=bool)
        d1[:] = tokeep
      #self.task = "PCA"
      #self.subtask = "FILTERS"


# PCA stuff 

  def PCA_filters(self, filter_list, value_list):
    if PCA_LD in filter_list: 
      ind = filter_list.index(PCA_LD)
      winsize = value_list.pop(ind)
      filter_list.pop(ind)
      self.r1 = int(winsize)
      self.r2 = int(winsize/2)
      n = 0
      for chrom in self.store:
        n = max(n, len(self.store["{}/positions".format(chrom)]))
      self.tqdm = tqdm.tqdm(total=n)

    self.QC_filters(filter_list, value_list, prefix="PCA_")


  def ld_filters(self, message, thresh=0.2):
    self.counter -= 1
    msg = {"TASK":"PCA", "SUBTASK":"LD"}
    for key, val in message.items():
      if key == "TASK" or key == "SUBTASK":
        continue
      else:
        if key in self.sumLin:
          self.sumLin[key] += val[0]
          self.sumSq[key]  += val[1]
          self.cross[key]  += val[2]
        else:
          self.sumLin[key] = val[0]
          self.sumSq[key]  = val[1]
          self.cross[key]  = val[2]
 
    if self.counter == 0: # Every report is in
      for key, val in message.items():
        if key == "TASK" or key == "SUBTASK":
          continue
        corr_tot = corr(self.sumLin[key], self.sumSq[key], self.cross[key])
        group  = self.store[key]
        tokeep = self.store["{}/PCA_passed".format(key)].value
        end = min(self.r1 + self.r0, len(tokeep))
        maf = self.store["{}/allele_freq".format(key)].value[self.r0:end]
        maf = maf[tokeep[self.r0:end]]
        n = maf.shape[0]
        unfiltered = np.ones((n, ), dtype=bool)
        while True:
          #length_of_window = np.sum(tokeep[self.r0:end)
          for i, snp1 in enumerate(unfiltered):
            if not snp1: # already filtered 
              continue
            else: 
              for j in range(i+1,n):# range(loc+1-self.r0, self.r1):
                snp2 = unfiltered[j]
                if not snp2: # if it didn't pass the filters
                  continue
                elif corr_tot[i,j] ** 2 > thresh:
                  if maf[i] > maf[j] * (1.0 + kLargeEpsilon):
                    unfiltered[i] = False
                  else: 
                    unfiltered[j] = False
                  break
          r2 = corr_tot[unfiltered,:][:,unfiltered]
          if np.max(r2**2) < thresh:
            break
        tokeep[self.r0:end][tokeep[self.r0:end]] = unfiltered
        pca_passed = self.store["{}/PCA_passed".format(key)]
        pca_passed[:] = tokeep
        end = min(self.r1 + self.r0 + self.r2, len(tokeep))
        if self.r1 + self.r0 >= len(tokeep):# 
          msg[key]="END"
        else:
          msg[key] = tokeep[self.r0 + self.r2:end]
      self.counter = self.connections
      self.r0 += self.r2
      self.sumLin, self.sumSq, self.cross = dict(), dict(), dict()
      self.server.message(encode(msg))
      self.tqdm.update(self.r2)


  def report_content(self, dset_name, msg):
    chroms = [v for v in self.store.keys() if v != 'meta']
    original_msg = msg.copy()
    for chrom in chroms: 
      dset = self.store["{}/{}".format(chrom, dset_name)]
      msg[chrom] = dset.value
      self.server.message(encode(msg))
      time.sleep(1)
      msg = original_msg.copy()

  def buildCov(message):
    # store the chunks in the store. build on them 
    ch1 = msg["CH1"] 
    ch2 = msg["CH2"]
    group = self.store["meta"]
    if "{}_{}".format(ch1, ch2) in group:
      del group["{}_{}".format(ch1, ch2)]
    group.create_dataset("{}_{}".format(ch1, ch2), data = message["MAT"])




if __name__=='__main__':
  print("no commands here yet. Test using WTCCC_run.py")


