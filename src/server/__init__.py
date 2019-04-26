import os
import logging
import sys
sys.path.append(os.path.abspath('../lib'))

FORMAT = '[{levelname:8}] {asctime:25} [server] {message}'
err_fmt = f'[%(levelname)-5.5s] %(asctime)s %(pathname)-80.80s :: %(lineno)s => %(message)s'
logging.basicConfig(level=logging.INFO, format=err_fmt, style='%')
