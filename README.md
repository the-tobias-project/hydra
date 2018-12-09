# HYDRA
Decentralized GWAS

## Data prep:

We will use chromosome 21 and 22 from 1000 Genomes phase 3 for practice (n=2504). You can find the relevant dataset in smallData or 
download a copy using smallData/download.sh For the purposes of this tutorial, the data has been thinned to contain 100k snps. 



## Server Setup:

py3 server.py


## Client Setup:

py3 runner.py

(runner.py starts 3 data silos) each data silo can be started using `py3 client.py  data_plink_stem`

## Performing a GWAS

When all data silos are connected you are prompted to specify how you'd like to proceed. 

`Looks like everyone is here. Type the stage you want to proceed with? Options: HELP, INIT, QC, PCA, EXIT:`
`init` 

During this time the files are read and some preliminary statistics are computed.
On the server side, you should see:
  Now working on: Initialized
  50223 loci in chromosome 21
  49777 loci in chromosome 22
  Now working on: INIT POS
  Count statistics have been initialized!
  Count statistics have been initialized!
  Count statistics have been initialized!
  Count statistics have been initialized!
  Count statistics have been initialized!
  Now working on: INIT STATS

Once the initial setup is over, The user is prompted to specify QC filters:
`Indicate the filters and corresponding values(e.g. hwe 10e-5). Available filters are HWE, MAF, MPS(Missing Per sample), MPN(Missing per snp), snp(rsid) (MPS not implemented yet):`
`hwe 10e-5`

This step is relatively quick and you should see the following response from the server: 
  in chromosome 21, 2833 snps were deleted and 47390 snps remain
  in chromosome 22, 3103 snps were deleted and 46674 snps remain

The initial command prompt reappears. At this point, we can exit to further analyze the data, we can apply more filters or proceed to perform PCA: 
`Looks like everyone is here. Type the stage you want to proceed with? Options: HELP, INIT, QC, PCA, EXIT:`
`pca`
`Specify filters for PCA pre-processing: (OPTIONS: HWE, MAF, MPS, as before as well as LD. e.g. maf 0.1 LD 50:`
`maf 0.1 ld 50`
