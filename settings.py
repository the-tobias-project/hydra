"""
Classes in this file are standalone because we don't want to impose a false hierarchy
between two classes.  That is, inheritance may imply a hierarchy that isn't real.
"""


class Settings(object):
    kExactTestBias = 1.0339757656912846e-25
    kSmallEpsilon = 5.684341886080802e-14
    kLargeEpsilon: 1e-07
    SMALL_EPSILON: 5.684341886080802e-14
    local_scratch = "/local/scratch/armin/Hydra"  # This should be stored in a more ephemeral location


class Commands(object):
    HELP = "HELP"
    INIT = "INIT"
    QC = "QC"
    PCA = "PCA"
    ASSO = "ASSO"
    EXIT = "EXIT"
    all_commands = [HELP, INIT, QC, PCA, ASSO, EXIT]


class Options(object):
    HELP = Commands.HELP
    INIT = Commands.INIT
    QC = Commands.QC
    PCA = Commands.PCA
    ASSO = Commands.ASSO
    EXIT = Commands.EXIT
    HWE = "HWE"
    MAF = "MAF"
    MPS = "MPS"
    MPI = "MPI"
    SNP = "snp"
    LD = "LD"


class QCOptions(object):
    HWE = Options.HWE
    MAF = Options.MAF
    MPS = Options.MPS
    MPI = Options.MPI
    SNP = Options.SNP
    all_options = [HWE, MAF, MPS, MPI, SNP]


class PCAOptions(object):
    HWE = Options.HWE
    MAF = Options.MAF
    MPS = Options.MPS
    MPI = Options.MPI
    SNP = Options.SNP
    LD = Options.LD
    all_options = [HWE, MAF, MPS, MPI, SNP, LD]


class QCHWE(object):
    QC_HWE = Options.HWE
    QC_MAF = Options.MAF
    QC_MPS = Options.MPS
    QC_MPI = Options.MPI
    QC_snp = Options.SNP


class PCAHWE(object):
    PCA_HWE = Options.HWE
    PCA_MAF = Options.MAF
    PCA_MPS = Options.MPS
    PCA_MPI = Options.MPI
    PCA_snp = Options.SNP
    PCA_LD = Options.LD