# stdlib
import logging

# Third party lib
from flask import Blueprint, request
from flask import current_app as app

# Internal lib
from lib import networking

bp = Blueprint('root', __name__, url_prefix='/api')
#celery_client = Celery(broker=Settings.redis_uri, backend=Settings.redis_uri)
from worker.tasks import celery as celery_client

@bp.route('/init', methods=['POST'])
def init():
    logging.info('Got command to initialize')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_store',
            [app.config['client'], app.config["ENV"]],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/init/stats', methods=['POST'])
def init_stats():
    logging.info('Got command to store initialized stats')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_stats',
                            [request.data, app.config['client'], app.config["ENV"]],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/delayed', methods=['GET'])
def delayed():
    logging.info('called delayed celery entry point')
    client_name = app.config['client']['name']
    promise = celery_client.send_task('tasks.caller', [adder_fn, 1, 2],
                                      serializer='pickle',
                                      queue=client_name)
    # resolution = promise.wait()  # answer is in here, if the celery backend is defined.  This will block.
    return networking.create_response(200)


@bp.route('/after_delayed', methods=['GET'])
def after_delayed():
    logging.info('called after_delayed celery entry point')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.dependent', None, serializer='pickle', queue=client_name)
    return networking.create_response(200)


@bp.route('/qc', methods=['POST'])
def qc():
    logging.info('Got command for QC')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.init_qc',
                            [request.data, app.config['client'], app.config['ENV']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/pca/ld', methods=['POST'])
def ld_report():
    logging.info('Got command for LD')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.report_ld',
                            [request.data, app.config['client'], app.config['ENV']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/pca/pcapos', methods=['POST'])
def store_filtered():
    logging.info('Got results of filtered positions')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.store_filtered',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/pca/cov', methods=['POST'])
def communicate_cov():
    logging.info('Preparing to report covariances')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.report_cov',
                            [app.config['client'], app.config['ENV']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/pca/eig', methods=['POST'])
def pca_projection():
    logging.info('Computing projections')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.pca',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/asso/adjust', methods=['POST'])
def data_adjust():
    logging.info('Covariate update')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.adjust',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/asso/initialize', methods=['POST'])
def lr_init():
    logging.info('Initializing Regression')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.regression_init',
                            [app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/asso/estimate', methods=['POST'])
def lr_association():
    logging.info('Regression update')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.asso',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/asso/coef', methods=['POST'])
def compute_likelihoods():
    logging.info('Regression update')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.loglikelihood',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)


@bp.route('/asso/query', methods=['POST'])
def compute_cost():
    logging.info('Performing line search')
    client_name = app.config['client']['name']
    celery_client.send_task('tasks.lineSearch',
                            [request.data, app.config['client']],
                            serializer='pickle',
                            queue=client_name)
    return networking.create_response(200)



def adder_fn(a, b):
    return a + b
