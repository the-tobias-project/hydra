#!/usr/bin/env python3

# stdlib
import argparse
import logging
import os
import sys
import traceback
import json
from pickle import UnpicklingError

# Third party lib
from twisted.internet import reactor, protocol, threads
from twisted.internet.defer import DeferredQueue, inlineCallbacks

# Internal lib
import clientSideAnalysis
from src.lib.settings import Settings
from utils import decode

HOST = 'localhost'
CLIENT_NAME = 'NO-NAME'
plink = None
local_scratch = None


class MyClient(protocol.Protocol):
    def __init__(self):
        self.name = CLIENT_NAME  # should get passed in as argument to constructor...
        self.cache = None
        self.setup_logger()

    def setup_logger(self):
        # logFormatter = logging.Formatter("%(asctime)s [%(threadName)-50.50s] [%(levelname)-5.5s]  %(message)s")
        logFormatter = logging.Formatter(f"%(asctime)s [{self.name:10}] [%(levelname)-5.5s]  %(message)s")
        self.logger = logging.getLogger()
        fileHandler = logging.FileHandler("{0}/{1}.log".format(os.getcwd(), "HYDRA_{}_logger".format(self)))
        fileHandler.setFormatter(logFormatter)
        self.logger.addHandler(fileHandler)
        consoleHandler = logging.StreamHandler(sys.stdout) 
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.setLevel(logging.INFO)

    def connectionMade(self):
        self.logger.info("Connection made!")
        self.factory.clients.append(self)
        self.logger.info("Clients are {}".format(self.factory.clients))
        self.cdispatcher = clientSideAnalysis.ServerTalker(plink, 
            local_scratch, self)
        self.wait_for_and_process_next_message()

    def clientConnectionLost(self, reason):
        # TODO send warning
        return()

    def dataReceived(self, data):
        if self.cache is not None:
            data = self.cache + data
        try:
            message = decode(data)
            message_queue.put(message)
            self.cache = None
        except json.decoder.JSONDecodeError as e:
            self.cache = data
        except UnpicklingError as e:
            self.cache = data
        except Exception as e:
            logging.error(e)
            logging.error(traceback.format_exc())
            self.cache = data

    def message(self, data):
        reactor.callFromThread(self._message, data)

    def _message(self, data):
        self.transport.write(data)

    def wait_for_and_process_next_message(self):
        reactor.callLater(0, self._wait_for_and_process_next_message)

    @inlineCallbacks
    def _wait_for_and_process_next_message(self):
        message = yield message_queue.get()
        yield threads.deferToThread(self.cdispatcher.dispatch, message)
        self.wait_for_and_process_next_message()


class MyClientFactory(protocol.ClientFactory):
    protocol = MyClient


def main():
    parser = argparse.ArgumentParser(description='CWS client')
    parser.add_argument('plinkfile', type=str, help='The plinkfile to process')
    parser.add_argument('--local_scratch', type=str, help='Location used for scratch storage during computation')
    parser.add_argument('--port', type=int, help='Which port the clients should use for communication')
    parser.add_argument('--name', type=str, help='A name for this client, useful in logging')
    args = parser.parse_args()

    global plink
    global local_scratch
    plink = args.plinkfile

    if args.local_scratch is not None:
        local_scratch = args.local_scratch
    else:
        local_scratch = Settings.local_scratch

    if args.name is not None:
        global CLIENT_NAME
        CLIENT_NAME = args.name
    factory = MyClientFactory()
    if args.port is not None:
        PORT = args.port
    else:
        PORT = 9000
    reactor.connectTCP(HOST, PORT, factory)
    factory.clients = []
    reactor.run()


message_queue = DeferredQueue()

if __name__ == "__main__":
    main()

