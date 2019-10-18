# stdlib
import logging
import logging.config

# Third party libs
from celery import Celery
from celery.utils.nodenames import gethostname

# internal libs
from lib.logging_config import return_worker_config

def make_celery(app_name=__name__):
    app = Celery(app_name)
    app.config_from_object("worker.settings")
    return app

celery = make_celery('cws_queue')
worker_config = return_worker_config(gethostname() + '.log')
logging.config.dictConfig(worker_config)
