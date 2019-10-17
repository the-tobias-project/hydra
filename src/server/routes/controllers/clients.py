# stdlib
import logging

# internal lib
from lib.client_registry import Registry
from lib import networking


def lst_clients():
    """
    List all registered clients
    """

    registry = Registry.get_instance()
    msg = registry.list_clients()
    return networking.create_response(200, msg)


def add_client(client):
    """
    Add a client to the registry
    """
    registry = Registry.get_instance()
    added = registry.add_client(client)
    if not added:
        return networking.create_response(400, 'Client already registered')
    logging.info(f"""Added {client['name']}, listening host: {client['listen_host']}, \
external_host: {client['external_host']}""")
    return networking.create_response(200)


def remove_client(client_name):
    """
    Remove a client from the registry
    """
    if not isinstance(client_name, str):
        msg = 'client name must be a string'
        return networking.create_response(400, msg)

    registry = Registry.get_instance()
    registry.remove_client(client_name)
    logging.info(f"Removed {client_name}")
    return networking.create_response(200)


def report_status(client_name, status):
    logging.info(f'[{client_name}]: {status}')
    return networking.create_response(200)
