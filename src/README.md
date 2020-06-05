# Directory structure

├── client  
│   ├── __pycache__  
│   ├── lib  
│   │   └── __pycache__  
│   └── routes  
│       └── __pycache__  
├── lib  
│   └── __pycache__  
├── server  
│   ├── __pycache__  
│   ├── lib  
│   │   └── __pycache__  
│   └── routes  
│       ├── __pycache__  
│       ├── controllers  
│       │   └── __pycache__  
│       └── schemas  
└── worker  
    └── __pycache__  

## client

Contains code related to the spoke side of the software.  
lib: helper functions for cleaning up paths.  
routes/dispatcher.py: Dispatches commands from the hub to the the celery queue. 

## lib 

This folder contains files that are useful for both hub and spoke. such as optimizationAux. There is probably some room for clean up here.  
Of particular importance are networking.py which governs how the hub and spoke talk to each other and settings which sets parameter values for both hub and spokes. 

## server 

Contains code related to the hub side of the software
routes/schemas has the swagger schema file.  
routes/controllers Contains files that manage response to the inputs from browser and the spokes. clients.py manages the clients (presence, deletion etc.) and tasks.py manages computational tasks. 

## worker

Contains code related to the spoke. The client code sends the jobs to these workers. 

Files of the form `task_*.py` correspond to the code that does the "work" the rest of the files are to set up the workers and register the tasks (`tasks.py`)






