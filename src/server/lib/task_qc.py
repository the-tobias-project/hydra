# stdlib
import logging
import pickle

# third party lib
import requests

# internal lib
from lib.settings import Settings
from lib.client_registry import Registry

clients = Registry.get_instance().list_clients()


def start_qc_task():
    for client in clients:
        Registry.get_instance().set_client_state(client['name'], 'QC')
        requests.post(f'http://{client["external_host"]}/{client["port"]}/api/qc')
