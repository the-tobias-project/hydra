# rmin Pourshafeie

#TODO documentation
import json, bson
import _pickle as pickle
import numpy as np
import plinkio, json, tqdm, time, sys
import h5py, os, tqdm, itertools
from plinkio import plinkfile
import pdb

# In house packages
from corr import nancorr, corr, HweP
#from optimizationAux import * # To be used later

with open("GLOBALS.json", 'r') as fp:
  locals_dict = json.load(fp)
for key, val in locals_dict.items():
  exec(key + '=val')


def encode(message):
  return pickle.dumps(message)

def decode(message):
  return pickle.loads(message)

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
    self.store_path      = self._store_path()
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

  def _store_path(self):
    plinkName = self.store_name
    basename = os.path.basename(plinkName)
    write_to = self.scratch
    store_name = os.path.join(write_to, basename + '.h5py')
    return store_name


  def plinkToH5(self):
    """Gets plink prefix, produces an HDF file with the same prefix"""
    plinkName = self.store_name
    plink_file = plinkfile.open(plinkName)
    if not plink_file.one_locus_per_row():
      print( "This script requires that snps are rows and samples columns.")
      sys.exit(1)
    
    sample_list = plink_file.get_samples()
    locus_list = plink_file.get_loci()
    n_tot = len(sample_list)
    with h5py.File(self.store_path, 'w', libver='latest', swmr=True) as store: 
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
      genotypes = np.zeros(n_tot, dtype=np.float32)
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
        genotypes[:] = np.array(row, dtype=np.float32)
        pos = str(locus.bp_position)
        dset = current_group.require_dataset(pos, (n_tot,), dtype=np.float32)
        positions.append(pos)
        dset.attrs['rsid'] = locus.name
        dset.attrs['snp']  = locus.allele1
        dset.attrs['alt']  = locus.allele2
        #vals, counts = np.unique(genotypes, return_counts=True)
        #dset.attrs['counts'] = np.array([counts[i] for i in range(len(vals)) if vals[i] != 0], dtype=np.int32)
        dset.attrs['counts'] = np.array([np.sum(genotypes==1), np.sum(genotypes==2), np.sum(genotypes==3)] , dtype=np.int32)
        genotypes[genotypes == 3] = np.nan
        dset[:] = genotypes 
        
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
        chrom = self.store_message(message, "PCA_passed")
        self.chroms.remove(chrom)
        if not self.chroms:
          self.chroms = None
          self.load_chroms()
          self.standardize_data(scale=False, center=True)
          self.report_covariance(self.chroms , "PCA_passed")



      pass
    elif task == "Association":
      pass



  def standardize_data(self, scale=False, center=True):
    with h5py.File(self.store_path, 'a') as store: 
      for chrom in self.chroms:
        group = store[chrom]
        if center: 
          af = group["MAF"].value
        if scale:
          sd = np.sqrt(group["VAR"].value)
        pos  = group["positions"].value
        for i, snp in enumerate(pos):
          genotypes = group[str(snp)]
          if center:
            vals = genotypes.value
            vals[np.isnan(vals)] = 2*af[i]
            genotypes[:] = vals - 2*af[i]
          if scale: 
            genotypes[:] = genotypes[:] / sd[i]


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
      if "VAR" in message:
        vals = message["VAR"]
        task = "var"
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
  def report_covariance(self, chroms, mask_name):
    msg = {"TASK": "PCA", "SUBTASK":"COV"}
    with h5py.File(self.store_path, 'r') as store:
      n = store.attrs["n"]
      for i_ch1, ch1 in enumerate(chroms):
        group = store[ch1]
        pos   = group["positions"].value
        mask  = group[mask_name].value
        pos   = pos[mask]
        g1     = np.empty((len(pos), n))

        for i, snp1 in enumerate(pos):
          g1[i, :] = group[str(snp1)].value

        for i_ch2, ch2 in enumerate(chroms):
          if i_ch2 > i_ch1 :
            break
          group = store[ch2]
          pos   = group["positions"].value
          mask  = group[mask_name].value
          pos   = pos[mask]
          g2     = np.empty((n, len(pos)))
          
          for i, snp2 in enumerate(pos):
            g2[:, i] = group[str(snp2)].value

          msg["CH1"] = ch1
          msg["CH2"] = ch2
          msg["MAT"] = g1.dot(g2).astype(np.float32)
          print(ch1, ch2)


   
if __name__=='__main__':
  print("no commands here yet. Test using WTCCC_run.py")


