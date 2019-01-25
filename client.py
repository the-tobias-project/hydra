#!/usr/bin/env python3

# stdlib
import time
import sys, os
import logging
import pdb

# Third party lib
from twisted.internet import reactor, protocol, threads
from twisted.internet.defer import DeferredQueue, inlineCallbacks

# Internal lib
import clientSideAnalysis
from serverSideAnalysis import decode
from settings import Settings

HOST = 'localhost'


class MyClient(protocol.Protocol):
    def __init__(self):
        self.cache = None
        self.setup_logger()

    def setup_logger(self):
        logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
        self.logger = logging.getLogger()
        fileHandler = logging.FileHandler("{0}/{1}.log".format(os.getcwd(), "HYDRA_client_logger"))
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
        except Exception as e:
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

message_queue = DeferredQueue()
if __name__=="__main__":
    plink = sys.argv[1]
    if len(sys.argv) >= 3:
        local_scratch = sys.argv[2]
        if len(sys.argv) == 4:
            PORT = int(sys.argv[3])
        else:
            PORT = 9000
    else:
        local_scratch = Settings.local_scratch

    factory = MyClientFactory()
    reactor.connectTCP(HOST, PORT, factory)
    factory.clients = []
    reactor.run()
