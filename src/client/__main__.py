# stdlib
import argparse
import logging
import signal
import sys

# third party lib
from flask import Flask
import requests

# Internal lib
from client.lib import shared
from lib import settings

# Routes
from client.routes import tasks

app = Flask(__name__)

server_host = settings.ServerHTTP.external_host
server_port = settings.ServerHTTP.port
BASE_URL = f'http://{server_host}:{server_port}/api'


def parse_args():
    client_names = list(map(lambda x: x['name'], settings.ClientHTTP.clients))

    parser = argparse.ArgumentParser(description='CWS client')
    parser.add_argument('--name',
                        type=str,
                        help='Name of the client to start.',
                        choices=client_names,
                        required=True)
    """
    The plinkfile is not of type argparse.FileType, as the base name won't exist
    
    e.g. popres1 is what we use to get popres1.[ind, bim, bam, fam], but popres1 doesn't exist.
    """
    parser.add_argument('--plinkfile',
                        type=str,
                        help='The plinkfile to analyze',
                        required=True)
    parser.add_argument('--port', type=int, help='[OPTIONAL] Override the default port')
    parser.add_argument('--external_host', type=str, help='[OPTIONAL] Override the default host used by '
                                                          'external systems to access this client.  Defaults to '
                                                          f'{settings.ClientHTTP.default_external_host}')
    parser.add_argument('--max_len', type=int, help='[OPTIONAL] Maximum content length for a given request.'
                                                    f'Defaults to {settings.ClientHTTP.default_max_content_length} b')
    parser.add_argument('--listen_host', type=str, help='[OPTIONAL] Override the default host on which this client'
                                                        'should listen.  Defaults to '
                                                        f'{settings.ClientHTTP.default_listen_host}')

    return parser.parse_args()


def setup_logging(client_name):
    # Sorry in advance - the below mixes two different mini formats!
    fmt = f'[%(levelname)-5.5s] %(asctime)s [{client_name:10}] %(message)s'
    logging.basicConfig(level=logging.INFO, format=fmt, style='%')
    err_fmt = f'[%(levelname)-5.5s] %(asctime)s [{client_name:10}] %(pathname)-100.100s :: %(lineno)s => %(message)s'
    logging.basicConfig(level=logging.ERROR, format=err_fmt, style='%')


def configure_client(client, args):
    if args.external_host is not None:
        client['external_host'] = args.external_host
    if args.port is not None:
        client['port'] = args.port
    if args.listen_host is not None:
        client['listen_host'] = args.listen_host
    if args.max_len is not None:
        client['max_content_length'] = args.max_len
    client['plinkfile'] = args.plinkfile
    shared.set_plinkfile(args.plinkfile)

    return client


def register_self(client):
    url = f'{BASE_URL}/clients'
    try:
        registered_clients = requests.get(url).json()
        self_name = client['name']
    except Exception as e:
        logging.error('Error getting list of registered clients from server')
        logging.error(e)
        return False
    try:
        # Already registered, no need to register again
        next(filter(lambda x: x['name'] == self_name, registered_clients['msg']))
        logging.info('Already registered with server, not attempting to register again.')
    except StopIteration:
        # Not registered
        requests.post(url, json=client)
        logging.info('Successfully registered self with server')
    except Exception as e:
        logging.error(e)
        return False
    return True


def teardown(signum, frame):
    try:
        client = app.config['client']
        url = f'{BASE_URL}/clients'
        requests.delete(f'{url}/{client["name"]}')
        sys.exit(0)
    except Exception as e:
        logging.error('Ran into unexpected error during teardown')
        logging.error(e)
        sys.exit(1)


def main():
    args = parse_args()
    client = next(filter(lambda x: x['name'] == args.name, settings.ClientHTTP.clients))
    setup_logging(args.name)

    app.config.update(
        CELERY_BROKER_URL='redis://localhost:6379',
        CELERY_RESULT_BACKEND='redis://localhost:6379'
    )

    # Overrides
    client = configure_client(client, args)

    app.config['MAX_CONTENT_LENGTH'] = client['max_content_length']
    app.config['client'] = client  # Store configuration for later use

    # Handle a teardown on sigint/sigterm
    signal.signal(signal.SIGINT, teardown)
    signal.signal(signal.SIGTERM, teardown)

    app.register_blueprint(tasks.bp)

    if not register_self(client):
        logging.error('Could not register self with server, exiting...')
        sys.exit(1)

    app.run(host=client['listen_host'], port=client['port'], threaded=False)


if __name__ == '__main__':
    main()
