import os
import sys

# stdlib
import socket
import time

# third party lib
#from flask_celeryext import create_celery_app
#from celery import current_app

# internal lib
from worker import task_init, task_qc, task_pca, task_asso

sys.path.append(os.path.abspath('../lib'))
sys.path.append(os.path.abspath('../client'))
#app = Flask('__name__')
#app.config["CELERY_BROKER_URL"] = Settings.redis_uri
#app.config["CELERY_BACKEND_URL"] = Settings.redis_uri
#app = create_celery_app(app)
#celery = Celery('cws_queue', broker=Settings.redis_uri, backend=Settings.redis_uri)
from worker import celery
from flask import current_app as app

@celery.task(name='tasks.hello')
def hello():
    return 'hello world'


@celery.task(name='tasks.caller')
def caller(fn, a, b):
    print(f'calling supplied function with two values {a} and {b}')
    result = fn(a, b)
    time.sleep(20)
    print(f'And the result is {result}')
    return result


@celery.task(name='tasks.dependent')
def dependent():
    """
    {'celery@ubuntu-bionic':
    [{'id': '814be1f1-be83-4f25-8331-cb9a1b7d4337', 'name': 'tasks.caller',
    'args': '[<function adder_fn at 0x7fd3247d2048>, 1, 2]', 'kwargs': '{}', 'type': 'tasks.caller',
    'hostname': 'celery@ubuntu-bionic', 'time_start': 1550479401.2852194, 'acknowledged': True,
    'delivery_info': {'exchange': '', 'routing_key': 'celery', 'priority': 0, 'redelivered': None},
    'worker_pid': 6974},
    {'id': '3a0c23e5-955a-4fc6-8983-2d1760d7dedf', 'name': 'tasks.dependent',
    'args': '()', 'kwargs': '{}', 'type': 'tasks.dependent', 'hostname': 'celery@ubuntu-bionic',
    'time_start': 1550479402.5393896, 'acknowledged': True, 'delivery_info': {'exchange': '',
    'routing_key': 'celery', 'priority': 0, 'redelivered': None}, 'worker_pid': 6973}]}
    """
    print(f'Called a dependent function')
    hostname = socket.gethostname()
    i = app.control.inspect()
    times_called = 0
    while True:
        times_called += 1
        active_tasks = i.active()[f'celery@{hostname}']
        dependent_tasks = list(filter(lambda x: x['type'] == 'tasks.caller', active_tasks))
        if times_called == 1:
            print('Remaining tasks that are still active:')
            print(dependent_tasks)
        if len(dependent_tasks) > 0:
            print('Waiting on tasks to finish...')
            time.sleep(1)
        else:
            break
    print('Broke free!')


@celery.task(name='tasks.init_store')
def init_store(client_config, env):
    task_init.init_store(client_config, env)


@celery.task(name='tasks.init_stats')
def init_stats(message, client_config, env):
    task_init.init_stats(message, client_config, env)


@celery.task(name='tasks.init_qc')
def init_qc(message, client_config, env):
    task_qc.init_qc(message, client_config, env)


@celery.task(name='tasks.report_ld')
def report_ld(message, client_config, env):
    ld_agg = task_pca.LdReporter.get_instance(50, client_config)
    ld_agg.update(message, client_config, env)


@celery.task(name='tasks.store_filtered')
def store_filtered(message, client_config):
    task_pca.store_filtered(message, client_config)


@celery.task(name='tasks.report_cov')
def report_cov(client_config, env):
    task_pca.report_cov(client_config, env)


@celery.task(name='tasks.pca')
def pca_projection(data, client_config):
    task_pca.pca_projection(data, client_config)


@celery.task(name='tasks.regression_init')
def initialize_logistic_reg(client_config, env):
    lr_agg = task_asso.LogisticAdmm.get_instance(range(2, 4), 10, client_config, env)


@celery.task(name='tasks.asso')
def compute_logistic_reg(message, client_config):
    print("Running regression")
    lr_agg = task_asso.LogisticAdmm.get_instance(range(2, 4), 10, client_config)
    lr_agg.update(message, client_config)


@celery.task(name='tasks.adjust')
def adjust_covariates(message, client_config, env):
    print("adjusting covariates")
    lr_agg = task_asso.LogisticAdmm.get_instance(range(2, 4), 10, client_config, env)
    lr_agg.global_standardize(message, client_config)


@celery.task(name='tasks.loglikelihood')
def compute_log_likelihood(message, client_config):
    lr_agg = task_asso.LogisticAdmm.get_instance(range(2, 4), 10, client_config)
    lr_agg.send_likelihood(message)


@celery.task(name='tasks.lineSearch')
def compute_cost(message, client_config):
    lr_agg = task_asso.LogisticAdmm.get_instance(range(2, 4), 10, client_config)
    lr_agg.cost(message)
