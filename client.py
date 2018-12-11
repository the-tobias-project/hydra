from twisted.internet import reactor, protocol, threads
import time
import clientSideAnalysis
from serverSideAnalysis import decode 
import sys
from twisted.internet.defer import DeferredQueue, inlineCallbacks
import pdb


HOST = 'localhost'
PORT = 9000
local_scratch="/local/scratch/armin/Hydra"

class MyClient(protocol.Protocol):
  def __init__(self):
    self.cache = None

  def connectionMade(self):
    print("connected!")
    self.factory.clients.append(self)
    print ("clients are ", self.factory.clients)

    self.cdispatcher = clientSideAnalysis.ServerTalker(plink, local_scratch, self)
    self.wait_for_and_process_next_message()


  def clientConnectionLost(self, reason):
    #TODO send warning
    self.factory.clients.remove(self)

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

#factory = MyClientFactory()
#reactor.connectTCP(HOST, PORT, factory)
#factory.clients = []
#reactor.run()
#
message_queue = DeferredQueue()
if __name__=="__main__":
  plink = sys.argv[1]

  factory = MyClientFactory()
  reactor.connectTCP(HOST, PORT, factory)
  factory.clients = []
  reactor.run()

