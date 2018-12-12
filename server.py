from twisted.internet import reactor, protocol, threads
from twisted.internet.defer import DeferredQueue, inlineCallbacks
import time
from serverSideAnalysis import ServerTalker, decode
import pdb, sys
import json
import _pickle as pickle


### SPLIT THE CONNECTION SO THAT IT CAN GO BACK TO OTHER TASKS
#TODO: graceful exit 

PORT = 9000
NUM = 3
local_scratch="/local/scratch/armin/Hydra"

with open("GLOBALS.json", 'r') as fp:
  local_dict = json.load(fp)
for key, val in local_dict.items():
  exec(key + '=val')
  
HELP_STRING="HELP MEEEEE"

class Hub(protocol.Protocol):
    def __init__(self,factory, clients, nclients):
      self.clients = clients 
      self.nclients = nclients
      self.factory = factory
      self.dispatcher = ServerTalker(NUM, self, 
          local_scratch)
      self.cache   = None

    def get_response(self, options):
      options = [val.upper() for val in options]
      waiting_for_command = True
      while waiting_for_command:
        val = input("Looks like everyone is here. Type the stage you want to proceed with? Options: " + ", ".join(options) 
            +": " )
        val = val.upper()
        if val in options:
          waiting_for_command = False
      if val == EXIT:
        self.tearDown()
      elif val == HELP:
        print(HELP_STRING)
      elif val == INIT: 
        self.run_init()
      elif val == "QC":
        self.run_QC()
      elif val == "PCA":
        self.run_pca()
        print(val)
      elif val == "ASS": 
        waiting_for_command = False
      self.wait_for_and_process_next_message()
    
    def run_QC(self):
      task = "QC"
      while True:
        val = input("""Indicate the filters and corresponding values(e.g. hwe 1e-5). Available filters are HWE, MAF, MPS(Missing Per sample), MPN(Missing per snp), snp(rsid) (MPS not implemented yet): """)
        val = val.upper()
        vals = val.split()
        if vals[0] == EXIT:
          self.tearDown()
          return
        elif len(vals) < 2:
          print("Specify `filter value` pairs? OK? OOOK?")
        else: 
          subtasks = [v.upper() for v in vals[::2]]
          vals     = vals[1::2]
          # Sanity checks:
          if not len(set(subtasks)) == len(subtasks):
            print("You can only use the same filter once.")
            continue
          if not set(subtasks).issubset(QC_OPTIONS):
            print("Wrong keywords")
            continue 
          break
      vals = [float(v) for v in vals]
      outMessage = pickle.dumps({"TASK": task, "SUBTASK": "FILTERING", "FILTERS":subtasks, "VALS":vals})
      self.message(outMessage)
      inMessage= { "TASK":"QC", "SUBTASK":subtasks, "VALS":vals}
      message_queue.put(inMessage)
      message_queue.put("GET OPTIONS")
      #reactor.callLater(2, self.run_pca())
      

    def run_init(self):
      message = pickle.dumps({"TASK": "INIT", "SUBTASK":"STORE"})
      self.message(message)
    
    def run_pca(self): 
      task = "PCA" 
      while True: 
        val = input("""Specify filters for PCA pre-processing: (OPTIONS: HWE, MAF, MPS, as before as well as LD. e.g. maf 0.1 LD 50: """)
        vals = val.split()
        if val == EXIT:
          sys.exit()
        if len(vals) < 2:
          print("Specify `filter value` pairs. PLEASEEEE!! REEEE!")
        else:
          subtasks = [v.upper() for v in vals[::2]]
          vals     = vals[1::2]
          assert len(set(subtasks)) == len(subtasks)
          assert set(subtasks).issubset(PCA_OPTIONS)
          break 
      vals = [float(v) for v in vals]
      outMessage = pickle.dumps({"TASK": task, "SUBTASK": "FILTERING", "FILTERS":subtasks, "VALS":vals})
      self.message(outMessage)
      inMessage = {"TASK":"PCA", "SUBTASK":"FILTERS", "FILTERS":subtasks, "VALS": vals}
      message_queue.put(inMessage)

    def connectionMade(self):
      print("connected to user" , (self))
      if len(self.clients) < self.nclients:
        self.factory.clients.append(self)
      else:
        self.factory.clients[self.nclients] = self
      if len(self.clients) == NUM:
          val = self.get_response(OPTIONS)


    def wait_for_and_process_next_message(self):
      reactor.callLater(0, self._wait_for_and_process_next_message)


    @inlineCallbacks
    def _wait_for_and_process_next_message(self):
#      if not message_queue.pending and self.moveOn:
#        self.get_response(OPTIONS)

      message = yield message_queue.get()
      if message == "GET OPTIONS":
        self.get_response(OPTIONS)
      else:
        yield threads.deferToThread(self.dispatcher.dispatch, message)
        self.wait_for_and_process_next_message()

    def message(self, command):
      reactor.callFromThread(self._message, command)

    def _message(self, command):
        for c in self.factory.clients:
          c.transport.write(command)

    def tearDown(self): 
      reactor.stop()
      for c in self.factory.clients:
        c.connectionLost("Exiting")
      return

    def connectionLost(self, reason): #TODO
      self.factory.clients.remove(self)
      self.nclients -= 1

    def dataReceived(self, data):
      if len(self.clients) == NUM:
        if self.cache is not None:
          data = self.cache + data
        try : 
          message = decode(data)
          message_queue.put(message)
          self.cache = None
        except Exception as e:
          self.cache = data


class PauseTransport(protocol.Protocol):
  def makeConnection(self, transport):
      transport.pauseProducing()

class HubFactory(protocol.Factory):
    def __init__(self, num):
      self.clients = set([])
      self.nclients = 0 
      self.totConnections = num

    def buildProtocol(self, addr):
      print(self.nclients)
      if self.nclients < self.totConnections:
        self.nclients += 1
        return Hub(self, self.clients, self.nclients)
      protocol = PauseTransport()
      protocol.factory = self
      return protocol

factory = HubFactory(NUM)
reactor.listenTCP(PORT, factory)
factory.clients = []
message_queue = DeferredQueue()
reactor.run()
