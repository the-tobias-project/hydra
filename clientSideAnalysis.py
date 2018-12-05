# rmin Pourshafeie

#TODO write a generator that takes the chromosome and spits out data. do the regression in parallel 
#TODO documentation
# Running the gwas
import logging
import pdb
import numpy as np
import gzip, re, gc
from sklearn.linear_model import LogisticRegression
import statsmodels.formula.api as smf
from statsmodels.tools.tools import add_constant
from functools import partial
from pathos.multiprocessing import ProcessingPool as Pool
import sklearn.decomposition as decomp 
from scipy.linalg import svd
from scipy.stats import chi2 
from scipy.sparse.linalg import eigsh as eig 
from optimizationAux import *
from plinkio import plinkfile
# Careful here, eigh uses https://software.intel.com/en-us/mkl-developer-reference-c-syevr behind the hood
# so it can be significantly slower 
from numpy.core import _methods
from sklearn.utils.extmath import randomized_svd, svd_flip
import time, sys
from numpy.linalg import inv as inverse
from numpy.core import umath as um
from numpy import mean, isnan
from sklearn.metrics import log_loss



import json, bson
import _pickle as pickle
import plinkio
import h5py, os, tqdm, itertools
from corr import nancorr, corr, HweP

with open("GLOBALS.json", 'r') as fp:
  locals_dict = json.load(fp)
for key, val in locals_dict.items():
  exec(key + '=val')


def encode(message):
  return pickle.dumps(message)

def decode(message):
  return pickle.loads(message)

#from numpy.core import umath as um
#umr_maximum = um.maximum.reduce
umr_sum = um.add.reduce
maximum = np.maximum
add     = np.add
_mean   = _methods._mean
_sum    = _methods._sum
sub     = np.subtract 
div     = np.divide 
chi2sf  = chi2.sf
sqrt    = np.sqrt

mean = np.mean

kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;
kSmallEpsilon = 0.00000000000005684341886080801486968994140625;
kLargeEpsilon = 1e-7

class ServerTalker(object):
  """Dispatches commands from/to the server"""
  def __init__(self, store_name, scratch, server):
    self.store_name = store_name
    self.status     = ""
    self.scratch    = scratch
    self.server     = server
    self.store_path      = None
    if not os.path.exists(scratch):
      os.makedirs(scratch)
    self.local_n = None

    self.chroms = None
    self.do_ld = False
    self.r1     = None
    self.r2     = None
    self.r3, self.r4 = None, None
    self.tqdm   = None
    self.verbose = True



  def plinkToH5(self):
    """Gets plink prefix, produces an HDF file with the same prefix"""
    plinkName = self.store_name
    plink_file = plinkfile.open(plinkName)
    write_to = self.scratch
    if not plink_file.one_locus_per_row():
      print( "This script requires that snps are rows and samples columns.")
      sys.exit(1)
    
    sample_list = plink_file.get_samples()
    locus_list = plink_file.get_loci()
    n_tot = len(sample_list)
    basename = os.path.basename(plinkName)
    store_name = os.path.join(write_to, basename + '.h5py')
    self.store_path = store_name
    with h5py.File(store_name, 'w', libver='latest', swmr=True) as store: 
      store.attrs['n'] = len(sample_list)
      store.attrs['has_local_AF'] = False
      store.attrs['has_global_AF'] = False
      store.attrs['has_centering'] = False
      store.attrs['has_normalization'] = False
      dset = store.require_dataset('meta/Status', (n_tot,), dtype=np.int8)
      dset[:,] = [sample_list[i].affection for i in range(n_tot)]
      # Ids
      dset = store.require_dataset('meta/id', (n_tot,), dtype='S11')
      dset[:,] = [sample_list[i].iid.encode('utf8')  for i in range(n_tot)]
      # Read Demographic file
      dem_file = plinkName + ".ind"
      with open(dem_file, 'r') as dem_f: 
        dset = store.require_dataset('meta/regions', (n_tot,), dtype='S19')
        dset[:,] = ["DatRegion".encode('utf8') for i in range(n_tot)] #[line.split('\t')[5] for line in dem_f]
      # Read chromosome data
      current_chr = 1
      positions = []
      current_group = store.require_group(str(current_chr))
      genotypes = np.zeros(n_tot, dtype=np.int8)
      for locus in locus_list:
        row = next(plink_file)
        if locus.chromosome != current_chr: 
          if len(positions) == 0:
            del store[str(current_chr)]
          else:
            dset = current_group.require_dataset('positions',
                (len(positions),), dtype=np.uint32) 
            positions = np.array(positions, dtype=np.uint32)
            dset[:] = np.array(positions, dtype=np.uint32)
            msg = encode({"TASK":"INIT", "SUBTASK":"POS", "CHROM":current_chr, "POS":positions})
            self.server.message(msg)
            positions = []
          current_chr = locus.chromosome
          if current_chr == 23:
            break
          current_group = store.require_group(str(current_chr))
        genotypes[:] = np.array(row, dtype=np.int8)
        pos = str(locus.bp_position)
        dset = current_group.require_dataset(pos, (n_tot,), dtype=np.float32)
        positions.append(pos)
        dset[:,] = genotypes 
        dset.attrs['rsid'] = locus.name
        dset.attrs['snp']  = locus.allele1
        dset.attrs['alt']  = locus.allele2
        #vals, counts = np.unique(genotypes, return_counts=True)
        #dset.attrs['counts'] = np.array([counts[i] for i in range(len(vals)) if vals[i] != 0], dtype=np.int32)
        dset.attrs['counts'] = np.array([np.sum(genotypes==1), np.sum(genotypes==2), np.sum(genotypes==3)] , dtype=np.int32)
      if locus.chromosome != 23:
        dset = current_group.require_dataset('positions',
            (len(positions),), dtype=np.uint32) 
        dset[:] = np.array(positions, dtype=np.uint32)
        msg = encode({"TASK":"INIT", "SUBTASK":"POS", "CHROM":current_chr, "POS":positions})
        self.server.message(msg)



  def dispatch(self,message):
    message = message
    task = message["TASK"]
    subtask = message["SUBTASK"]
    if self.verbose:
      print("Working on:", task)
    self.status = "{}/{}".format(task, subtask) 
    if task == "INIT":
      if subtask == "STORE":
        self.plinkToH5()
        print("preparing counts")
        self.reportCounts()
      if subtask == "STATS":
        self.store_stats(message)
    elif task == "QC":
      self.run_QC(message)
    elif task == "PCA":
      # First we will perform all the filters save for LD pruning: 
      if "FILTERS" in message: 
        self.run_snp_filters(message)
      elif self.do_ld: # Just do LD
        self.run_ld(message)
      elif subtask == "PCA_POS":
        if self.chroms is None:
          self.load_chroms()
        #chrom = [key for key, _ in message.items() if key not in ["TASK", "SUBTASK"]][0]
        #if str(chrom) in self.chroms: 
          #self.chroms.remove(str(chrom))
        self.store_message(message, "PCA_passed")
        if not self.chroms:
          self.chroms = None
          self.load_chroms()
          self.center_data()
          self.report_covariance(self.chroms , "PCA_passed")



      pass
    elif task == "Association":
      pass

  def center_data(self):
    with h5py.File(self.store_path, 'a') as store: 
      for chrom in self.chroms:
        af   = group["{}/MAF".format(chrom)]
        dset = group[chrom]
        pos  = group["{}/positions".format(chrom)]
        
      



  def load_chroms(self):
    with h5py.File(self.store_path, 'r') as store:
      self.chroms = [i for i in store.keys() if i!= 'meta']


  def store_stats(self, message):
    chrom = message["CHROM"]
    with h5py.File(self.store_path, 'a') as store: 
      chrom_group = store[str(chrom)]
      if "MISS" in message: 
        vals = message["MISS"]
        task = "Missing_per_snp"
        dset = chrom_group.create_dataset(task, data=vals)
      if "AF" in message:
        vals = message["AF"]
        task = 'MAF'
        dset = chrom_group.create_dataset(task, data=vals)
      if "HWE" in message:
        vals = message["HWE"]
        task = "hwe"
        dset = chrom_group.create_dataset(task, data=vals)


  def reportCounts(self):
    """Report the counts (Het, homo Alt, missing)"""

    with h5py.File(self.store_path, 'r') as store:
      countDict = {}
      n = store.attrs["n"]
      countDict["START"] = True
      keys = [i for i in store.keys() if i != 'meta']
      for chrom in keys:
        countDict["n"] = int(n)
        countDict["TASK"] ="INIT"
        countDict["SUBTASK"] = "COUNT"
        countDict["CHROM"] = chrom
        if chrom == 'meta':
          continue
        positions = store["{}/positions".format(chrom)].value
        count_arr = []
        for i, pos in enumerate(positions):
          count_arr.append(store["{}/{}".format(chrom, pos)].attrs["counts"])
        countDict["COUNTS"] = np.array(count_arr, dtype=np.uint32)#.tolist()
        if chrom == keys[-1]:
          countDict["END"] = True

        msg = encode(countDict)
        self.server.message(msg)
        countDict = {}



####### QC
  def run_QC(self, message, remove=True):
    def find_what_passes(qc_name, dset_name, tokeep, doubleSided=False):
      vals = group[dset_name].value
      if qc_name in filter_list:
        ind = filter_list.index(qc_name)
        if not doubleSided:
          tokeep = np.logical_and(tokeep, vals > value_list[ind])
        else:
          tokeep = np.logical_and(tokeep, np.logical_and(vals > value_list[ind], (1.0-vals) > value_list[ind]))
      return tokeep

    def replace_dataset(tokeep, dset_name, return_deleted=False):
      vals = group[dset_name].value
      del group[dset_name]
      remaining = vals[tokeep]
      deleted   = vals[np.logical_not(tokeep)]
      group.create_dataset(dset_name, data=remaining)
      if return_deleted:
        return deleted


    filter_list = message["FILTERS"]
    value_list   = message["VALS"]
    with h5py.File(self.store_path, 'a') as store:
      if self.chroms is None:
        self.chroms = [i for i in store.keys() if i != "meta"]

      for chrom in self.chroms:
        group = store[chrom]
        positions = group['positions'].value
        tokeep    = np.ones_like(positions, dtype=bool)
        tokeep    = find_what_passes(QC_HWE, "hwe", tokeep)
        tokeep    = find_what_passes(QC_MAF, "MAF", tokeep, doubleSided=True)
        tokeep    = find_what_passes(QC_MPS, "Missing_per_snp", tokeep)
        print("After filtering, {} snps remain".format(np.sum(tokeep)))
        if remove: # Delete what doesn't pass 
          replace_dataset(tokeep, 'hwe')
          replace_dataset(tokeep, 'MAF')
          replace_dataset(tokeep, 'Missing_per_snp')
          deleted = replace_dataset(tokeep, 'positions',return_deleted=True)
          for snp in deleted:
            del group[str(snp)]
        else: # Store what has been tagged
          if "passed" in group: 
            tags = group["passed"]
            tokeep = np.logical_and(tokeep, tags.value)
            tags[:] = tokeep
          else:
            group.create_dataset("passed", data=tokeep)

  ### PCA
  def run_snp_filters(self, message):
    print("Preparing for PCA")
    filters = message["FILTERS"]
    vals    = message["VALS"]
    if PCA_LD in filters: 
      ind = filters.index(PCA_LD)
      winsize = vals.pop(ind)
      filters.pop(ind)
      self.do_ld = True # deal with LD after other filters
      self.r1    = int(winsize)
      self.r2    = int(winsize/2)
    self.run_QC({"FILTERS":filters, "VALS":vals}, remove=False)
    if self.do_ld:
      self.r3 = 0
      with h5py.File(self.store_path, 'r') as fp: 
        if self.chroms is None:
          self.chroms = [i for i in store.keys() if i!= 'meta']
        n = 0
        for chrom in self.chroms:
          n = max(n, len(fp["{}/positions".format(chrom)]))

      print("LD pruning. This will take some time...")
      self.tqdm = tqdm.tqdm(total=n)
      self.verbose = False
      self.run_ld({})



  #def run_ld(self, winsize, step_size):
  def run_ld(self, message):
    msg = {"TASK": "PCA", "SUBTASK":"LD"}
    with h5py.File(self.store_path, 'r') as store:
      n = store.attrs['n']
      for chrom in self.chroms:
        tags = store["{}/passed".format(chrom)]
        if chrom in message:
          state = message[chrom]
          if state[0] == "E": # Finished with this chrom
            self.chroms.remove(chrom)
            if len(self.chroms) == 0:
              self.do_ld = False 
              self.chroms = None
              msg = {"TASK":"PCA", "SUBTASK":"PCA_POS"}
              self.server.message(encode(msg))
              self.verbose = True
              return
            continue
          else: 
            tokeep = state
            end = self.r3 + len(tokeep)
        else:
          end = min(self.r3 + self.r1, tags.shape[0])
          tokeep = tags[self.r3: end]
 
        pos = store["{}/positions".format(chrom)]
        positions = pos[self.r3: end]
        positions = positions[tokeep] 
        genotypes = np.empty((n,len(positions)), dtype=np.float32)
        for i, snp in enumerate(positions): 
          genotypes[:,i] = store["{}/{}".format(chrom, snp)].value
        corr = nancorr(genotypes)
        msg[chrom] = corr
    self.server.message(encode(msg))
    self.r3 += self.r2
    self.tqdm.update(self.r2)

  def store_message(self, message, dset_name):
    for key, val in message.items():
      if key in ["TASK", "SUBTASK"]:
        continue 
      else: 
        with h5py.File(self.store_path, 'a') as store:
          dset = store.require_dataset("{}/{}".format(key, dset_name), val.shape, dtype = val.dtype)
          dset[:] = val
#          corr = self.compute_corr(key, val)
#          msg = {"TASK":"PCA", "SUBTASK":"CORR", "CORR": corr, "LOC":(key, key)}
#          self.server.message(encode(msg))
        return key
  def report_covariance(self, chroms, mask):
    msg = {"TASK": self.task, "SUBTASK":"COV"}
    with h5py.File(self.store_path, 'r') as store:
      for ch1, ch2 in itertools.product(chroms, chroms): 
        group1 = store[ch1]
        pos1   = group1["positions"].value
        mask1  = group1[mask].value
        group2 = store[ch2]
        pos2   = group2["positions"].value
        mask2  = group2[mask].value








   

class Center(object):
  """The central hub that drives and requires particular computations from each node."""
  def __init__(self, server, store_names, n_cores=1):
    self.store_names = server.clients
    self.nDOs = len(store_names)
    self.ncores = n_cores
    # get preliminary Statistics
    package = {"Task": "Chr Report"}
    server.communicate(msgpack.packb(package))
    self.keys = self.DOs[0].group_keys()
    self.n   = sum([item.n for item in self.DOs])
    logging.info("- Setup center with {} DOs and {} individuals". format(
      self.nDOs, self.n))

  def loci_missing_rate_filter(self, rate):
    for DO in self.DOs: 
      DO.compute_local_missing_rates()
    for chrom in self.keys:
      if chrom != 'meta':
        logging.info("Consensus missing rate computation on chrom: {}".format(chrom))
        MR = add_dict()
        MR.set_key_values(self.DOs[0].dataset_keys(chrom), 0)
        for DO in self.DOs:
          update_dic = DO.report_local_missing_rates(chrom)
          MR.update(update_dic, 1.0)
        for DO in self.DOs: 
          DO.set_missing_rate_filter(chrom, MR, rate * self.n)


  def MAF_filter(self, rate):
    """Computes local and consensus AF, sd. 
    Removes loci below the specified MAF"""
    def AF_wrapper(DO):
      DO.compute_local_AF()
    #with Pool(self.ncores) as pool:
      #pool.map(AF_wrapper , self.DOs)
    for DO in self.DOs:
      AF_wrapper(DO)
    for chrom in self.keys:
      if chrom != 'meta':
        logging.info("---Consensus AF computation on chrom: {}".format(chrom))
        AF = add_dict()
        AF.set_key_values(self.DOs[0].dataset_keys(chrom),[0,0,0])
        for DO in self.DOs:
          update_dic = DO.report_local_AF(chrom)
          AF.update(update_dic, 1.0, 0) 
        # update the overall AF
        for DO in self.DOs: 
          DO.set_local_AF(chrom, AF, 0)
          if rate is not None:
            DO.MAF_filter(chrom, rate)

  def compute_std(self, chrom):
    if chrom != 'meta':
      logging.info("--consensus SD computation on chrom: {}".format(chrom))
      SD = add_dict()
      SD.set_key_values(self.DOs[0].dataset_keys(chrom), [0,0])
      for DO in self.DOs:
        update_dic = DO.report_SD(chrom)
        SD.update(update_dic, 1.0, 1) #TODO This is a colossal fuck up (AF,SD, HWE). All of this shit needs to be done by passing counts As the sufficient statistics. but too late for that shit now. will clean up later
      for DO in self.DOs: 
        DO.set_local_AF(chrom, SD, 1)


  def normalize(self):
    for chrom in self.keys:
      if chrom != 'meta':
        logging.info("--normalizing chrom: {}".format(chrom))
        self.compute_std(chrom)
        for DO in self.DOs: 
          DO.normalize(chrom)


  def HWE_filter(self, rate):
    for chrom in self.keys:
      if chrom != 'meta':
        logging.info("-HWE computation on chrom: {}".format(chrom))
        HWE = add_dict()
        HWE.set_key_values(self.DOs[0].dataset_keys(chrom),np.array([0,0,0]))
        for DO in self.DOs:
          update_dic = DO.report_local_counts(chrom)
          HWE.update(update_dic, 1.0)
        for key, value in HWE.iteritems(): 
          hwe = HweP(int(value[1]), int(value[0]), int(value[2]), 0 )
          HWE[key] = hwe
        for DO in self.DOs: 
          DO.HWE_filter(chrom, HWE, rate)
  
  def HWE_test(self, homor, het, homoa):
    """HWE test (midpoint test). Other versions of HWE filter can be impelemented with the same information. 
    This implementation should match PLINK1.9's implementation."""
    homc = max(homor, homoa)
    homr = min(homor, homoa)
    
    rare = 2 * homr + het
    # mid point of the distribution 
    n = (homor + het + homoa) * 2
    tail_p = (1 - kSmallEpsilon) * kExactTestBias
    centerp = 0
    lastp2, lastp1 = tailp, tailp 
    #if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)):

    mid = int(rare * (2 *  n -rare) / (2 * n))
    if (mid % 2 != rare % 2):
      mid += 1
    probs = np.zeros(1 + rare)
    probs[mid] = 1.0
    tsum = 1.0
    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = n - curr_hets - curr_homr
    
    while (curr_hets >= 2):
      probs[curr_hets - 2] = probs[curr_hets ] * (curr_hets) * (curr_hets - 1.0) / (4.0 * (curr_homr - 1.0) * (curr_homc + 1.0))
      tsum += probs[curr_hets - 2]
      curr_hets -= 2 
      curr_homr += 1
      curr_homc += 1
      
    curr_hets = mid 
    curr_homr = (rare - mid) / 2
    curr_homc = n - curr_hets - curr_homr
    while (curr_hets <= rare -2):
      probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      tsum += probs[curr_hets + 2]
      curr_hets += 2
      curr_homr -= 1
      curr_homc -= 1
    
#    target = probs[het]
#    return min(1.0, np.sum(probs[probs <= target])/tsum)

    probs /= tsum
    p_hi = float(probs[het])
    for i in xrange(het + 1, rare + 1):
      p_hi += probs[i]
#    
    p_lo = float(probs[het])
    for i in xrange(het-1, -1, -1):
      p_lo += probs[i]
    p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo
    p_hwe = 0.0
    for i in xrange(0, rare + 1):
      if probs[i] > probs[het]:
        continue
      p_hwe += probs[i]
    p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

    return p_hwe


  def correct_LD_prune(self, threshold, win_sz, step_sz=None):
    #TODO use local_LD_filter 
    def pruner(chrom, threshold, window):
      window.shape = (1, window.shape[0])
      to_delete = set()
      n = window.shape[1]
      sumLinT = np.zeros((n,n), dtype = np.float32)
      sumSqT = np.zeros((n,n), dtype = np.float32)
      crossT = np.zeros((n,n), dtype = np.float32)
      for DO in self.DOs: 
        sumLin, sumSq, cross = DO.corr_data([chrom], window)
        sumLinT += sumLin
        sumSqT  += sumSq
        crossT  += cross
      MAF = DO.get_MAF(chrom, window[0], global_freq=True)
      corrT = corr(sumLinT, sumSqT, crossT)
      while (1):
        for i, snp1 in enumerate(window[0,:]):
          if snp1 in to_delete:
            continue
          else:
            for j in range(i+1, n):
              if window[0][j] in to_delete:
                continue
              elif corrT[i,j]**2 > threshold:
                if MAF[i] > MAF[j] * (1.0 + kLargeEpsilon):   #somewhat similar to what plink does
                #ai = sumLin[i,j] / cross[i, j]
                #aj = sumLin[j,i] / cross[i, j]
                #majori = ai if ai > .5 else 1 - ai
                #majorj = aj if aj > .5 else 1 - aj
                #if ai > aj * (1 + kSmallEpsilon):
                  to_delete.add(snp1)
                else: 
                  to_delete.add(window[0][j])
                break
        remaining = np.array([i for i,snp in enumerate(window[0]) if snp not in to_delete])
        r2 = corrT[remaining,:][:,remaining]
        if np.max(r2**2) < threshold:
          break

      return to_delete

    if step_sz is None:
      step_sz = int(win_sz/2)
    for chrom in self.keys:
      if chrom == 'meta':
        continue 
      logging.debug("---Decentralized LD pruning on chrom: {}".format(chrom))
      # Get snps that pass the allele frequency threshold
      snps = np.sort(np.array(self.DOs[0].snps_present(chrom)).astype(int))
      win_sz = min(snps.shape[0], win_sz)
      finished, winstart  = False, 0
      highLD, to_delete = set(), set()
      while not finished:
        winend = winstart + win_sz
        if winend >= len(snps):
          finished = True 
          winend = len(snps)
        window = snps[winstart:winend] #preliminary window 
        window = np.sort(np.array(list(set(window) - to_delete)))#[:win_sz]
        to_delete = pruner(chrom, threshold, window)
        highLD = highLD.union(to_delete)
        winstart += step_sz# + offset[0][0]
      #toKeep = set(snps) - highLD
      logging.info("---Keeping {} snps after AF/LD pruning".format(len(snps) - len(highLD)))
      for DO in self.DOs:
        DO.delete_snps(chrom, highLD)

   

  def LD_prune(self,threshold, AF_threshold, win_sz, step_sz=None):
    """Flag snps that have small LD"""
    
    def pruner(chrom, threshold, window):
      window.shape = (1, window.shape[0])
      to_delete = set()
      n = window.shape[1]
      cov = np.zeros((n,n))
      # considerable optimization can be done so that only the parts 
      # that are previously not communicated get communicated 
      for DO in self.DOs:
        cov += float(DO.n)/float(self.n) * DO.give_cov([chrom], window)
      #cov /= self.nDOs
      # with covariance matrix we can be more accurate than the 
      # simple greedy we implemented in centralized but we go with the
      # same algorithm for comparison's sake 
      for i, snp in enumerate(window[0,:]):
        if snp in to_delete:
          continue
        else:
          for j in range(i+1, window.shape[1]):
            if window[0,j] in to_delete:
              continue
            elif cov[i,j]**2 > threshold:
              to_delete.add(window[0,j])
      return to_delete

    if step_sz == None:
      step_sz = int(win_sz/2)
    for chrom in self.keys:
      if chrom == 'meta':
        continue 
      logging.info("---Decentralized LD pruning on chrom: {}".format(chrom))
      # Get snps that pass the allele frequency threshold
      snps = np.sort(np.array(self.DOs[0].AF_filter(AF_threshold, chrom))).astype(int)
      win_sz = min(snps.shape[0], win_sz)
      finished, winstart  = False, 0
      highLD = set()
      i = 0
      while not finished:
        winend = winstart + win_sz
        if winend >= len(snps) - 1:
          finished = True 
          winend = len(snps) - 1
        window = snps[winstart:winend]
        window = np.sort(np.array(list(set(window) - highLD)))
        to_delete = pruner(chrom, threshold, window)
        highLD = highLD.union(to_delete)
        winstart += step_sz
        if winstart / 5000 > i:
          logging.debug("pruning at {}".format(winstart))
          i += 1
      toKeep = set(snps) - highLD
      logging.info("----Keeping {} snps after AF/LD pruning".format(len(toKeep)))
      for DO in self.DOs:
        DO.tag_snps(chrom, toKeep, 'prune_selected', True)

  def PCA(self, n_components=None, chroms=None):
    if chroms is None or chroms == []:
      chroms = [item for item in self.keys if item != 'meta']
    chroms = sorted(chroms, key=lambda x: int(x))
    DO = self.DOs[0]
    n = DO.count(list(set(self.keys) - set(chroms)))
    to_PCA = np.zeros((n, n), dtype=np.float32)
    logging.info("Preparing covariance matrix of size {}".format(n))
    for DO in self.DOs:
      DO.give_cov_pca(chroms, n, to_PCA, 1.0)# float(DO.n)/float(DO.n-1))
    if n_components is not None:
      m = min(self.n, n)
      m = min(m, n_components)
      #n_components = (n - n_components, n-1)
    #sigma, v = eig(to_PCA, overwrite_a=True, eigvals=n_components)# for linalg.eigh slow
    logging.info("Running PCA")
    sigma, v = eig(to_PCA, k=n_components, ncv=3*n_components)
    logging.info("Done running PCA")
    # there should be no ev with negativ e ev. If there is it should 
    # be tiny and due to numerical errors 

    del to_PCA
    sigma, v = zip(*sorted(zip(sigma, v.T),reverse=True))
    v = np.array(v)

    sigma = np.array(sigma)
    sigma[sigma < 0] = 0
    for DO in self.DOs:
      DO.store_eigs(sigma, v, chroms)
    #pca = PCA(n_components=n_components)
     #for now ignore the n_components arg
    #pca.fit(to_PCA)

  def change_pheno(self, pheno_plink):
    pheno_file = plinkfile.open(pheno_plink)
    sample_list = pheno_file.get_samples()
    iid = [item.iid for item in sample_list]
    status = [item.affection  for item  in sample_list]
    status_dict = dict((key, value) for (key, value) in zip(iid, status))
    for DO in self.DOs: 
      DO.update_pheno(status_dict)

  def copy_pca(self, other, local=False):
    for DO in self.DOs:
      base = os.path.basename(DO.store_name)
      file_name = os.path.join(other, base)
      DO.copy_pca(file_name, local)

  def run_regression(self, numPCs, n_iters, warm_start=True, chroms=[], sites=None, kind='ADMM',
      verbose=False, out_file="d_beta.txt"):

    def _regression(kind, verbose, **kwargs):
      """Dispatches to regression algorithm"""
      if kind == 'ADMM':
        if verbose:
          return self._ADMM_verbose(**kwargs)
        else:
          return self._ADMM(**kwargs)
      elif kind == 'AVG':
         return self._AVG(**kwargs)


    logging.info("-Running regression")
    DOs = self.DOs

    kwargs = {"rho": 10.0, "max_iters":n_iters, "alpha":1.2,
      "npcs":numPCs, "mu":0.0}#self.n * 1e-9}
    # Compute the variance of PCs
    first_moment = np.zeros((1, numPCs))
    second_moment = np.zeros((1, numPCs))
    #covp = len(pos) + numPCs
    covp = numPCs + 1

    for DO in DOs:
      DO.load_snp = True
      m1, m2 = DO.give_moments(["meta/pca_u"])
      first_moment += np.array(m1[0][:numPCs]) * DO.n / float(self.n)
      second_moment += np.array(m2[0][:numPCs]) * DO.n / float(self.n)
    stds = np.sqrt(second_moment - first_moment**2)
    kwargs["stds"] = stds

    write_n = 50
    if verbose:
      write_n = write_n / 10

    # Run for each snp
    if len(chroms) == 0 :
      chroms = self.keys
    else:
      chroms = [unicode(str(chrom)) for chrom in chroms]
    num_g = DOs[0].count(exclude=list(set(self.keys) - set(chroms)))
    pbar = tqdm.tqdm(total=num_g)
    counter, all_betas, warm_beta = 0, [], np.zeros((covp, 1))

    # Run regression with PC's only one time, to get the likelihood for the smaller model
    kwargs['pos'] = []
    kwargs["beta"] = warm_beta[1:]
    pc_beta = _regression(kind, False, **kwargs)
    pc_likelihood = 0
    
    warm_beta[1:] = pc_beta

    for DO in DOs: 
      pc_likelihood += DO.likelihood(pc_beta)
      DO.load_snp = True
      DO.current_Y = None
    
    if not verbose:
      pval = np.empty((covp + 2, 1))
    else:
      pval = np.empty((covp + 2, n_iters+1))
 
    # Run regression for everything else and compute the log likelihood difference/Wald Pvalues
    with open(out_file, 'w') as fout:
      for chrom in chroms:
        if chrom == 'meta':
          continue
        logging.info("--Running {} on chromosome: {}".format(kind, chrom))
        snps = sorted(DOs[0].dataset_keys(chrom), key=lambda x:int(x))
        pval[covp+1, :] = chrom
        for snp in snps:
          kwargs["pos"] = [(chrom, snp)]
          kwargs["beta"] = warm_beta
          beta = _regression(kind, verbose, **kwargs)
          if isnan(beta[0,0]):
            pval[:covp+1,:] = np.nan
            for DO in DOs:
              DO.load_snp = True
          else:
            likelihood = 0
            for DO in DOs:
              likelihood += DO.likelihood(beta, verbose)
            covLogit = _sum([DO.covLogit([(chrom, snp)], beta, stds, True) for DO in DOs], axis=0)
            # get pvalues
            covLogit = inverse(covLogit)
            z = (beta / sqrt(np.diag(covLogit)).reshape(covp, 1))
            z = z * z
            pval[:covp,:] = chi2sf(z, 1)
            pval[covp,:] = likelihood - pc_likelihood
          if not verbose:
            all_betas.append( "\t".join(map(str, beta[:,0])) +"\t" + "\t".join(map(str, pval[:,0]))) 
          else:
            for ind, line in enumerate(beta.T):
              all_betas.append( "\t".join(map(str, line)) +"\t" + "\t".join(map(str, pval[:,ind].tolist() + [ind])))


          counter += 1
          if counter == write_n:
            fout.write('\n'.join(all_betas))
            fout.write('\n')
            counter = 0
            all_betas = []
            pbar.update(write_n)
      fout.write('\n'.join(all_betas))

  def _ADMM(self, pos, npcs, rho, beta, alpha=1., max_iters=10, mu=0.0, stds=1, #1e-9, stds = 1,
      logistic=True, verbose=True): # mu is really self.n * mu
    """Performs ADMM regression. So far, only logistic regression is implemented."""
    DOs = self.DOs
    covp = len(pos) + npcs
    K = len(DOs)
    z = np.zeros((covp, K))
    u = np.zeros((covp, K))
#    shrink_param = mu / float(rho * K)
    for k in xrange(max_iters):
      for i, DO in enumerate(DOs): # can be parallelized 
        try:
          # z update:
          z[:,i] = DO.admm_update(pos, npcs,u[:,i, None], beta, rho, z[:,i, None], stds, logistic, covp)
        except ValueError:
          beta *= np.nan
          return beta
      # Update betas
      z_hat = add(alpha * z, sub(1.0, alpha) * beta)
#      meanVal = div(_sum(add(z_hat, u), 1)[:,None], K)
#      beta = div(_sum(add(z_hat, u), 1)[:,None], K)
      beta = div(umr_sum(z_hat,1)[:,None], K)
#      beta = sub(maximum(0, sub(meanVal, shrink_param)), maximum(0, -add(meanVal, shrink_param)))
      # Update u:
      u += sub(z_hat, beta)
    return beta


  def _ADMM_verbose(self, pos, npcs, rho, beta, alpha=1.0, max_iters=10, mu=1e-9, stds=1, 
      logistic=True):
    """Same as _ADMM except records the beta after every iteration. _ADMM avoids checking the 
    condition over and over again. Probably a stupid optimization but w/e"""
    DOs = self.DOs
    covp = len(pos) + npcs
    K = len(DOs)
    z = np.zeros((covp, K))
    u = np.zeros((covp, K))
    shrink_param = mu / float(rho * K)
    Betas = np.empty((covp, max_iters+1))
    Betas[:,0] = 0
    Us = np.empty((1, max_iters+1))
    Us[0,0] = 0
    for k in xrange(max_iters):
      for i, DO in enumerate(DOs): # can be parallelized 
        try:
          # z update:
          z[:,i] = DO.admm_update(pos, npcs,u[:,i, None], beta, rho, z[:,i, None], stds, logistic, covp)
        except ValueError:
          Betas[k+1:, :] = np.nan
          return beta
      # Update betas
      z_hat = add(alpha * z, sub(1.0, alpha) * beta)
      #meanVal = div(_sum(add(z_hat, u), 1)[:,None], K)
      #beta = sub(maximum(0, sub(meanVal, shrink_param)), maximum(0, -add(meanVal, shrink_param)))
      beta = div(umr_sum(add(z_hat, u), 1)[:,None], K)
      Betas[:,k+1] = beta[:,0]
      # Update u:
      u += sub(z_hat, beta)
      Us[0,k+1] = np.linalg.norm(u)
    return Betas

  def _AVG(self, pos, npcs, stds = 1, logistic=True, verbose=True, **kwargs): 
    """Performs Average regression. So far, only logistic regression is implemented.
    Performs the regression on de-centralized data. This simply averages all the results,
    for the actual analysis, we used inverse variance weighted averaging FE model."""
    covp = len(pos) + npcs
    DOs = self.DOs
    N = float(self.n)

    beta = np.zeros((covp, 1))
    for i, DO in enumerate(DOs): # can be parallelized 
#      try:
        beta += DO.run_regression(pos, npcs, beta, stds, logistic, covp).T * DO.n / N
#      except ValueError:
#        beta *= np.nan
#        return beta
    # Update betas
    return beta



  def PCA_Centralized(self, n_components=None, chroms=None):
    from sklearn.decomposition import PCA 
    if chroms is None or chroms == []:
      chroms = [item for item in self.keys if item != 'meta']
    chroms = sorted(chroms, key=lambda x: int(x))
    DO = self.DOs[0]
    n = DO.count(list(set(self.keys) - set(chroms)))
    data = np.empty((self.n, n), dtype=np.float32)
    logging.info("centralizing data just to run centralized PCA")
    start = 0
    for DO  in self.DOs: 
      data[start:start+DO.n,:] = DO.give_data(chroms,n)
      start += DO.n

    pca = decomp.PCA()
    U, S, V = pca._fit_truncated(data, n_components=n_components, svd_solver = 'arpack')
#    u, sigma, vt = randomized_svd(data, n_components, transpose=False)
#    u,vt = svd_flip(u, vt, u_based_decision=False)
    self.DOs[0].record_centralized_pca(S, U) 
    logging.info("Done with centralized PCA")

  def run_meta_filters(self, t_missing=None, t_AF=None, t_hwe=None, t_LD=None, win_sz=50, global_clean=False):
    def count(global_clean):
      unfiltered = 0
      for chrom in self.keys:
        if chrom == 'meta':
          continue
        present = self.DOs[0].locally_unfiltered(chrom)
        for DO in self.DOs[1:]:
          present = present.intersection(DO.locally_unfiltered(chrom))
        unfiltered += len(present)
        if global_clean: 
          for DO in self.DOs:
            DO.clean_by_local_filter(chrom, present)
      return(unfiltered)
    if t_missing is not None:
      logging.info("Starting local missing filter")
      for DO in self.DOs: 
        DO.local_missing_filter(t_missing)
      unfiltered = count(global_clean)
      logging.info("After missing rate filter {} snps remain".format(unfiltered))
    if t_AF is not None:
      logging.info("Starting local AF")
      for DO in self.DOs:
        DO.local_AF_filter(t_AF)
      unfiltered = count(global_clean)    
      logging.info("After AF filter {} snps remain".format(unfiltered))
    if t_hwe is not None:
      logging.info("Starting HWE filter")
      for DO in self.DOs:
        DO.local_HWE_filter(t_hwe)
      unfiltered = count(global_clean)
      logging.info("After HWE filter {} snps remain".format(unfiltered))
    if t_LD is not None:
      logging.info("Running LD filter")
      for DO in self.DOs:
        DO.local_LD_filter(t_LD, win_sz) #implement
      unfiltered = count(global_clean)
      logging.info("After LD filter {} snps remain".format(unfiltered))

  def run_local_pca(self, n_components=10, chroms=None):
    for DO in self.DOs: 
      DO.local_pca(n_components, chroms)

  def run_meta_regression(self, numPCs, out_file):
    logging.info("Starting meta regression...")
    chroms = self.keys
    with open(out_file, 'a') as fout:
      for chrom in chroms:
        if chrom == 'meta': 
          continue
        logging.info("Moving on to chrom " + chrom)
        for i, DO in enumerate(self.DOs):
          betas, standard_errors, pvals = DO.local_regression(numPCs, chrom)
          if not i: # first DO
            to_write = np.empty((len(betas), 3*len(self.DOs)+1))
          to_write[:,i] = betas[:,0]
          to_write[:,i+len(self.DOs)] = standard_errors[:,0]
          to_write[:,i+2*len(self.DOs)] = pvals[:,0]
        to_write[:,3*len(self.DOs)] = chrom
        np.savetxt(fout, to_write)
    logging.info("Finished Meta-regressions")




  def impute (self):
    for DO in self.DOs: 
      DO.impute()
      logging.info("DUUUUDE")

class add_dict(dict):
  def set_key_values(self, keys=None, value=None):
    if keys is None:
      keys = self.keys()
    if value is None:
      value = 0
    for key in keys:
      self[key] = value


  def update(self, other, frac=1.0, pos=None):
    if pos is None:
      k1 = other.keys()[0]
      if isinstance(other[k1], int):
        for key, value in other.iteritems():
          dicVal = self[key]
          self[key] = dicVal + frac * value
      else:# it is an array
        for key, value in other.iteritems():
          dicVal = self[key]
          self[key] = [x + frac * y for x,y in zip(dicVal, value)]
    elif pos == 0: #deal with these later TODO they can be put in the framework above
      for key, value in other.iteritems():
        dicVal = self[key]
        self[key] = dicVal[0] + value[2] * value[0], dicVal[1], dicVal[2] + value[2]
    elif pos == 1:
      for key, value in other.iteritems():
        dicVal = self[key]
        self[key] = dicVal[0] + value[0]**2, dicVal[1] + value[1]

if __name__=='__main__':
  print("no commands here yet. Test using WTCCC_run.py")


