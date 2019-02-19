# stdlib
import logging

# third party lib
from flask import request

# internal lib
from lib import tasks
from lib import HTTPResponse
from lib.settings import Commands
from server.lib import task_init
from server.lib import task_qc


def list_tasks():
    tsks = tasks.task_list
    return HTTPResponse.create_response(200, tsks)


def start_task(task_name):
    logging.info(f'Got command to start {task_name}, starting...')
    if task_name == Commands.INIT:
        task_init.start_init_task()
    elif task_name == Commands.QC:
        task_qc.start_qc_task()


def start_subtask(task_name, subtask_name, client_name):
    logging.info(f'Got task {task_name}/{subtask_name}')
    if task_name == Commands.INIT:
        if subtask_name == 'POS':
            logging.info(f'Got POS response from {client_name}')

            task_init.store_positions(request.data)
        elif subtask_name == 'COUNT':
            logging.info('Got a count subtask')
            logging.info(f'Got COUNT response from {client_name}')
            task_init.store_counts(request.data, client_name)

    return HTTPResponse.create_response(200)