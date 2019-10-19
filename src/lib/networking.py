# third party lib
from flask import jsonify, make_response
import requests

# internal lib
from lib.settings import ServerHTTP
from lib.client_registry import Registry


def create_response(code, msg=None):
    """
    Convenience method for flask.make_response(flask.jsonify(msg), code)

    If no message is supplied, a generic message will be returned based on the code provided

    :param code:            The HTTP status code
    :param msg:             Any message to be included in the response
    """
    status = 'OK'
    if code != 200:
        status = 'error'
    if msg is None:
        if code == 200:
            msg = 'success'
        elif code == 404:
            msg = 'Not found'
        elif code == 400:
            msg = 'User error; please check request'
        else:
            msg = 'Unknown error'
    envelope = {
        'status': status,
        'msg': msg
    }
    return make_response(jsonify(envelope), code)


def respond_to_server(path, verb, msg=None, client_name=None, env='production'):
    url = f'{get_protocol(env)}://{ServerHTTP.external_host}:{ServerHTTP.port}/{path}'
    s = requests.Session()
    req = requests.Request(method=verb, url=url, data=msg, params={'client_name': client_name})
    prepped = req.prepare()
    s.send(prepped)


def get_protocol(env='production'):
    if env == 'development':
        return 'http'
    return 'https'


def message_clients(address, client_name=None, args=None, env='production', data=None):
    clients = Registry.get_instance().list_clients()
    if client_name is None:
        client_list = clients
    else:
        client_list = list(filter(lambda x: x["name"] == client_name, clients))
    for client in client_list:
        if args is None:
            requests.post(f'{get_protocol(env)}://{client["external_host"]}:{client["port"]}/api/{address}',
                          data=data)
        else:
            requests.post(f'{get_protocol(env)}://{client["external_host"]}:{client["port"]}/api/{address}',
                          params=args, data=data)
