import os

from flask import Flask

# internal libs
from worker.celery_utils import init_celery

default_name = os.path.dirname(os.path.realpath(__file__)).split("/")[-1]
def create_app(name=default_name, **kwargs):
    application = Flask(name)
    if kwargs.get("celery"):
        init_celery(kwargs.get("celery"), application)
    from client.routes.tasks import bp
    application.register_blueprint(bp)
    return application


