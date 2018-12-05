from twisted.internet import reactor, protocol, threads
import msgpack
import time
import clientSideAnalysis
from serverSideAnalysis import decode 
import sys
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

  def clientConnectionLost(self, reason):
    #TODO send warning
    self.factory.clients.remove(self)

  def dataReceived(self, data):
    if self.cache is not None: 
      data = self.cache + data
    try:
      message = decode(data)
      self.cache = None
      threads.deferToThread(self.cdispatcher.dispatch, message)
    except Exception as e:
      self.cache = data

  def message(self, data):
    reactor.callFromThread(self._message, data)

  def _message(self, data):
    self.transport.write(data)

class MyClientFactory(protocol.ClientFactory):
    protocol = MyClient

#factory = MyClientFactory()
#reactor.connectTCP(HOST, PORT, factory)
#factory.clients = []
#reactor.run()
#
#message_queue = DeferredQueue()
if __name__=="__main__":
  plink = sys.argv[1]

  factory = MyClientFactory()
  reactor.connectTCP(HOST, PORT, factory)
  factory.clients = []
  reactor.run()

