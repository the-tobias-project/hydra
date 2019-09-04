"""
Classes in this file are standalone because we don't want to impose a false hierarchy
between two classes.  That is, inheritance may imply a hierarchy that isn't real.
"""


class Settings(object):
    kExactTestBias = 1.0339757656912846e-25
    kSmallEpsilon = 5.684341886080802e-14
    kLargeEpsilon = 1e-07
    SMALL_EPSILON = 5.684341886080802e-14
    local_scratch = '/app/scratch'
    python = 'python'
    plink = "/srv/gsfs0/software/plink/1.90/plink"
    redis_uri = 'redis://hydra_redis:6379'


class Commands(object):
    HELP = "HELP"
    INIT = "INIT"
    INIT_STATS = 'INIT_STATS'
    QC = "QC"
    PCA = "PCA"
    ASSO = "ASSO"
    EXIT = "EXIT"
    all_commands = [HELP, INIT, QC, PCA, ASSO, EXIT]  # used by v.1 interface


class Options(object):
    # HELP = Commands.HELP
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
    NONE = "NONE"


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
    NONE = Options.NONE
    all_options = [HWE, MAF, MPS, MPI, SNP, LD, NONE]


class QCFilterNames(object):
    QC_HWE = Options.HWE
    QC_MAF = Options.MAF
    QC_MPS = Options.MPS
    QC_MPI = Options.MPI
    QC_snp = Options.SNP


class PCAFilterNames(object):
    PCA_HWE = Options.HWE
    PCA_MAF = Options.MAF
    PCA_MPS = Options.MPS
    PCA_MPI = Options.MPI
    PCA_snp = Options.SNP
    PCA_LD = Options.LD
    PCA_NONE = Options.NONE


class ServerHTTP(object):
    listen_host = '0.0.0.0'
    external_host = 'localhost'
    port = '9001'
    max_content_length = 1024 * 1024 * 1024  # 1 GB
    wait_time = 0.5  # for the time.sleep() hacks


class ClientHTTP(object):
    default_max_content_length = 1024 * 1024 * 1024  # 1 GB
    default_listen_host = '0.0.0.0'
    default_external_host = 'localhost'
    clients = [{
            'name': 'BioME',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9002,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'MEC_CA',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9003,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'MEC_HI',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9004,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'SOL_B',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9005,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'SOL_C',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9006,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'SOL_M',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9007,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'SOL_S',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9008,
            'max_content_length': default_max_content_length
        },
        {
            'name': 'WHI',
            'listen_host': default_listen_host,
            'external_host': default_external_host,
            'port': 9009,
            'max_content_length': default_max_content_length
        }


    ]
