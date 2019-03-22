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
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/init/stats', methods=['POST'])
def init_stats():
    logging.info('Got command to store initialized stats')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_stats',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/delayed', methods=['GET'])
def delayed():
    logging.info('called delayed celery entry point')
    client_name = app.config['client']['name']
    promise = celery_client.send_task('tasks.caller', [adder_fn, 1, 2],
                                      serializer='pickle',
                                      queue=client_name)
    # resolution = promise.wait()  # answer is in here, if the celery backend is defined.  This will block.
    return HTTPResponse.create_response(200)


@bp.route('/after_delayed', methods=['GET'])
def after_delayed():
    logging.info('called after_delayed celery entry point')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.dependent', None, serializer='pickle', queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/qc', methods=['POST'])
def qc():
    logging.info('Got command for QC')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_qc',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/pca/ld', methods=['POST'])
def ld_report():
    logging.info('Got command for LD')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.report_ld',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/pca/pcapos', methods=['POST'])
def store_filtered():
    logging.info('Got results of filtered positions')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.store_filtered',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/pca/cov', methods=['POST'])
def communicate_cov():
    logging.info('Preparing to report covariances')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.report_cov',
                            [app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


@bp.route('/pca/eig', methods=['POST'])
def pca_projection():
    logging.info('Computing projections')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.pca',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return HTTPResponse.create_response(200)


def adder_fn(a, b):
    return a + b
