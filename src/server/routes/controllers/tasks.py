# stdlib
import logging

# third party lib
from flask import request

# internal lib
from lib import tasks
from lib import HTTPResponse
from lib.settings import Commands
from server.lib import task_init, task_qc, task_pca
from lib.client_registry import Registry


def list_tasks():
    tsks = tasks.task_list
    return HTTPResponse.create_response(200, tsks)


def start_task(task_name):
    logging.info(f'Got command to start {task_name}, starting...')
    if task_name == Commands.INIT:
        task_init.start_init_task()
    elif task_name.startswith(Commands.QC):
        task_name = "QChwe1e-10"
        filters = task_qc.split_command(task_name)
        print("Specified Filters :{filters}")
        task_qc.start_client_qc_task(filters)
        task_qc.start_local_qc_task(filters)
    elif task_name.startswith(Commands.PCA):
        if not task_pca.decomposed():
            if not task_pca.filtered():
                task_name = "PCAMAF0.1LD50_0.2"
                filters = task_qc.split_command(task_name)
                print(f"Specified pruning filters :{filters}")
                task_pca.start_pca_filters(filters)
            else:
                for client in Registry.get_instance().list_clients():
                    start_subtask(Commands.PCA, "PCAPOS", client)
        pass
    elif task_name == Commands.ASSO:
        pass


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

    elif task_name.startswith(Commands.QC):
        if subtask_name == "FIN":
            if task_qc.filter_finished(client_name, Commands.QC):
                print("We can move on")

    elif task_name.startswith(Commands.PCA):
        if subtask_name == "FIN":
            if task_qc.filter_finished(client_name, Commands.PCA):
                print("Done with PCA filters")
                reset_states("PRUNE")
                ld_agg = task_pca.CovarianceAggregator.get_instance(len(Registry.get_instance().list_clients()), 50)
        elif subtask_name == "LD":
            ld_agg = task_pca.CovarianceAggregator.get_instance(len(Registry.get_instance().list_clients()), 50)
            ld_agg.update(request.data)
        elif subtask_name == "PCAPOS":
            task_pca.report_pos(client_name)
            message_clients("pca/cov")


    return HTTPResponse.create_response(200)


def reset_states(state):
    instance = Registry.get_instance()
    for client in instance.list_clients():
        instance.set_client_state(client["name"], state)


def message_clients(message):
    for client in clients:
        requests.post(f'http://{client["external_host"]}:{client["port"]}/api/{message}')



