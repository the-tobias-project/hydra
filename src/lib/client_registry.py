# stdlib
import logging


class Registry(object):
    """
    Singleton, keeps track of registered clients
    """
    __instance = None

    @staticmethod
    def get_instance():
        if Registry.__instance is None:
            Registry()
        return Registry.__instance

    def __init__(self):
        if Registry.__instance is not None:
            # Short circuit the instantiation
            return
        else:
            self.registered_clients = []
            Registry.__instance = self

    def list_clients(self):
        return self.registered_clients

    def add_client(self, client):
        try:
            # If client is already registered, don't register it.
            next(filter(lambda x: x['name'] == client['name'], self.registered_clients))
            return False
        except StopIteration:
            self.registered_clients.append(client)
            return True

    def remove_client(self, client_name):
        clients = list(filter(lambda x: x['name'] != client_name, self.registered_clients))
        self.registered_clients = clients

    def set_client_state(self, client_name, state):
        try:
            client = next(filter(lambda x: x['name'] == client_name, self.registered_clients))
            client['state'] = state
        except StopIteration:
            logging.error(f'Cannot set client state for {client_name} - client not registered.')

    def get_client(self, client_name):
        try:
            return next(filter(lambda x: x['name'] == client_name, self.registered_clients))
        except StopIteration:
            logging.error(f'Cannot get client {client_name} - client not registered')
            return None

    def num_clients_in_state(self, state):
        return len(list(filter(lambda x: x['state'] == state, self.registered_clients)))
