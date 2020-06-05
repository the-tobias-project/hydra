# stdlib
import logging, pdb

# Third party lib


# Internal lib
from lib import networking
from worker.tasks import celery as celery_client
#from flask import current_app as app # want to eventually get rid of this

def dispatch_on_task(func):
    registry = {}

    def dispatch(value):
        try:
            return registry[value]
        except KeyError:
            return func

    def register(value, func=None):
        if func is None:
            return lambda f: register(value, f)
        registry[value] = func
        return func

    def wrapper(*args, **kw):
        return dispatch(args[0])(*args, **kw)

    wrapper.register = register
    wrapper.dispatch = dispatch
    wrapper.registry = registry

    return wrapper

@dispatch_on_task
def dispatcher(task, *args):
    pass
        


@dispatcher.register(("INIT", "start"))
def init(task, client, env, *args, **kw):
    logging.info('Got command to initialize')
    celery_client.send_task('tasks.init_store', [client, env],
            serializer='pickle', queue=client["name"])


@dispatcher.register(("INIT", "stats"))
def store_counts(task, client, env, *args, **kw):
    logging.info('Got command to store initialized stats')
    celery_client.send_task('tasks.init_stats',
                            [args[0], client, env],
                            serializer='pickle',
                            queue=client["name"])
