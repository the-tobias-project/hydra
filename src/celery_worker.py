from worker import celery
from worker.factory import create_app
from worker.celery_utils import init_celery

app = create_app()
init_celery(celery, app)
