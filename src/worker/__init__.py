
# Third party libs
from celery import Celery


def make_celery(app_name=__name__):
    app = Celery(app_name)
    app.config_from_object("worker.settings")
    return app


celery= make_celery('cws_queue')
