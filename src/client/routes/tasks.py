# stdlib
import logging

# Third party lib
from flask import Blueprint, request
from flask import current_app as app
from celery import Celery

# Internal lib
from lib import HTTPResponse
from lib.settings import Settings

bp = Blueprint('root', __name__, url_prefix='/api')
celery_client = Celery(broker=Settings.redis_uri, backend=Settings.redis_uri)


@bp.route('/init', methods=['POST'])
def init():
    logging.info('Got command to initialize')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_store',
                            [app.config['client']],
                            serializer='pickle',
                            headers={'client_name': client_name})
    return HTTPResponse.create_response(200)


@bp.route('/init/stats', methods=['POST'])
def init_stats():
    logging.info('Got command to store initialized stats')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_stats',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            headers={'client_name': client_name})
    return HTTPResponse.create_response(200)


@bp.route('/delayed', methods=['GET'])
def delayed():
    logging.info('called delayed celery entry point')
    promise = celery_client.send_task('tasks.caller', [adder_fn, 1, 2],
                                      serializer='pickle')
    # resolution = promise.wait()  # answer is in here, if the celery backend is defined.  This will block.
    return HTTPResponse.create_response(200)


@bp.route('/after_delayed', methods=['GET'])
def after_delayed():
    logging.info('called after_delayed celery entry point')
    celery_client.send_task('tasks.dependent', None, serializer='pickle')
    return HTTPResponse.create_response(200)


def adder_fn(a, b):
    return a + b