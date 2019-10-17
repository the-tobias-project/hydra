
# stdlib
import argparse
import io
import logging
import logging.config
import os
import yaml

# third party lib
import connexion
from flask import current_app

# Internal lib
from lib import settings
from lib.logging_config import return_config


options = {"swagger_ui": True}
app = connexion.FlaskApp(__name__, options=options)
BASE_SCHEMA = 'base.yml'

def parse_args():
    parser = argparse.ArgumentParser(description='CWS server')
    parser.add_argument('--port', type=int, help='[OPTIONAL] Override the default port')
    parser.add_argument('--external_host', type=str, help='[OPTIONAL] Override the default host used by '
                                                          'external systems to access this server.  Defaults to '
                                                          f'{settings.ServerHTTP.external_host}')
    parser.add_argument('--max_len', type=int, help='[OPTIONAL] Maximum content length for a given request.'
                                                    f'Defaults to {settings.ServerHTTP.max_content_length} b')
    parser.add_argument('--listen_host', type=str, help='[OPTIONAL] Override the default host on which this server'
                                                        'should listen.  Defaults to '
                                                        f'{settings.ServerHTTP.listen_host}')
    parser.add_argument('--scratch', type=str, help='[OPTIONAL] Override the default scratch location on disk.'
                                                    f'Defaults to {settings.Settings.local_scratch}')
    parser.add_argument('--dev', type=bool, default=False, help='[OPTIONAL] Specify a development environment.  '
                                                                'WARNING: this will bypass security checks.')
    parser.add_argument('--verbose', type=bool, default=False, help='[OPTIONAL] Specify log verbosity. ')
    return parser.parse_args()


def configure_server(server, args):
    if args.external_host is not None:
        server['external_host'] = args.external_host
    if args.port is not None:
        server['port'] = args.port
    if args.listen_host is not None:
        server['listen_host'] = args.listen_host
    if args.max_len is not None:
        server['max_content_length'] = args.max_len

    return server


def load_schemas():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    schema_dir = [dir_path, 'routes', 'schemas']
    schema_files = filter(lambda x: x.endswith('.yml'), os.listdir(os.path.join(dir_path, *schema_dir)))

    with io.open(os.path.join(*schema_dir, BASE_SCHEMA), 'r') as base_api:
        api = yaml.load(base_api, Loader=yaml.SafeLoader)

    for file in schema_files:
        if file == BASE_SCHEMA:
            continue

        with io.open(os.path.join(*schema_dir, file), 'r') as extension_file:
            api_ext = yaml.load(extension_file)

        if 'definitions' in api_ext:
            api['definitions'] = {**api['definitions'], **api_ext['definitions']}
        if 'paths' in api_ext:
            api['paths'] = {**api['paths'], **api_ext['paths']}
    try:
        outfile = 'generated.api.yml'
        with io.open(os.path.join(dir_path, outfile), 'w') as generated_api_file:
            yaml.dump(api, generated_api_file)
            return outfile
    except Exception as e:
        current_app.schema_info['status'] = 'error'
        current_app.schema_info['message'] = str(e)
        return False


def main():
    args = parse_args()
    logging_config = return_config(args.verbose)
    logging.config.dictConfig(logging_config)
    server = {
        'listen_host': settings.ServerHTTP.listen_host,
        'external_host': settings.ServerHTTP.external_host,
        'port': settings.ServerHTTP.port,
        'max_content_length': settings.ServerHTTP.max_content_length
    }
    server = configure_server(server, args)

    api_file = load_schemas()
    if api_file:
        app.add_api(api_file)
    else:
        logging.error('Error in initializing swagger specification')

    app.app.config['MAX_CONTENT_LENGTH'] = server['max_content_length']
    app.app.config['server'] = server  # Store configuration for later use
    if args.dev:
        app.app.config['ENV'] = 'development'

    app.run(host=server['listen_host'], port=server['port'], threaded=False)


if __name__ == '__main__':
    main()
