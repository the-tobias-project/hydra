# HYDRA
Decentralized GWAS

## Data prep:

We will use chromosome 21 and 22 from 1000 Genomes phase 3 for demonstration (n=2504). You can find the relevant dataset in testData or 
download a copy using testData/download1kG.sh For the purposes of this demo, the data has been thinned to contain 100k snps.And the individuals are split into 3 separate datasets. 



## Server Setup:
Out of the box, the server requires no arguments.  Assuming you're in the `src`
directory:
```bash
python3 -m server
```

## Infrastructure
We use [redis](https://redis.io/) to back the [celery](http://www.celeryproject.org/) workers - you
will need to have redis installed and running on the client machine(s).


## Client Setup:
The client needs a few additional arguments to launch successfully.  Assuming 
you're in the `src` directory, you can get a help message by running: 
```bash
python3 -m client --help
usage: __main__.py [-h] --name
                   {zealot,dragoon,high_templar,dark_templar,archon,dark_archon}
                   --plinkfile PLINKFILE [--port PORT]
                   [--external_host EXTERNAL_HOST] [--max_len MAX_LEN]
                   [--listen_host LISTEN_HOST]

CWS client

optional arguments:
  -h, --help            show this help message and exit
  --name {zealot,dragoon,high_templar,dark_templar,archon,dark_archon}
                        Name of the client to start.
  --plinkfile PLINKFILE
                        The plinkfile to analyze
  --port PORT           [OPTIONAL] Override the default port
  --external_host EXTERNAL_HOST
                        [OPTIONAL] Override the default host used by external
                        systems to access this client. Defaults to localhost
  --max_len MAX_LEN     [OPTIONAL] Maximum content length for a given
                        request.Defaults to 1073741824 b
  --listen_host LISTEN_HOST
                        [OPTIONAL] Override the default host on which this
                        clientshould listen. Defaults to 0.0.0.0
```
Although the `name` and `plinkfile` arguments are under the `optional arguments`
heading, attempting to run the client without these will print an error and throw
you back into your shell.

Let's say we have chosen to start the `zealot` client, whose plinkfiles are in
`/vagrant/testData/popres1`, i.e., there exist the following files:
 
 * `/vagrant/testData/popres1.bim`
 * `/vagrant/testData/popres1.bam`
 * (...)

Then to start the client, assuming the server is already running, from the `src` 
directory we run:

```bash
python3 -m client --name=zealot --plinkfile=/vagrant/testData/popres1
```

We also need a worker for the client, this can be set up by calling:

```bash
# Format: celery -A worker worker -Q <client_name>
# -A tells celery to look in the 'worker' directory
# The second 'worker' argument tells celery to start a worker

celery -A worker worker -Q zealot


  -------------- celery@ubuntu-bionic v4.2.1 (windowlicker)
---- **** -----
--- * ***  * -- Linux-4.15.0-45-generic-x86_64-with-Ubuntu-18.04-bionic 2019-02-20 19:58:42
-- * - **** ---
- ** ---------- [config]
- ** ---------- .> app:         cws_queue:0x7f12b2aa7d68
- ** ---------- .> transport:   redis://localhost:6379//
- ** ---------- .> results:     redis://localhost:6379/
- *** --- * --- .> concurrency: 2 (prefork)
-- ******* ---- .> task events: OFF (enable -E to monitor tasks in this worker)
--- ***** -----
 -------------- [queues]
                .> zealot           exchange=zealot(direct) key=zealot
```

You'll notice in Celery's initialization message above that it has registered itself to the
`zealot` queue.

If you look at the client and server logs, you'll notice the client registered 
itself with the server on initialization.  The server uses this registration to
later send out further tasks to each individual client.  On client SIGTERM or SIGINT,
(e.g. `control + c`) it will attempt to unregister itself from the server.

## Performing a GWAS
Assuming you want `N` clients, once you see that the server has `N` clients registered,
you can start the initialization step by calling:

```bash
curl http://localhost:9001/api/tasks/INIT -X POST
``` 

Where the server location and port number are set inside `src/lib/settings.py :: ServerHTTP`.

As the clients perform their tasks, they may send further tasks to the server - these 
will show up in the server logs, along with the name of the client that sent the task.
Additionally, once the task is completed, the client will send a status update to the
server, which will again show up in the server logs.  


### Server-side UI
If you're the kind of person that prefers a web-interface, first make sure you install
the UI-bindings for connexion by running:

```bash
pip3 install connexion[swagger-ui]
```

Then, assuming the server is on `localhost` on port `9001`, we can access the UI at:

[http://localhost:9001/api/ui/](http://localhost:9001/api/ui/)

From here you can click on `GET` requests to, say, get a list of registered clients,
or write in task paramters and then click a button to send a `POST` request to the
server.

### Deprecated?
When all data silos are connected you are prompted to specify how you'd like to proceed. 

`Looks like everyone is here. Type the stage you want to proceed with? Options: HELP, INIT, QC, PCA, EXIT:`  
`>init` 

During this time the files are read and some preliminary statistics are computed.

On the server side, you should see:  
```  
Now working on: Initialized  
50223 loci in chromosome 21  
49777 loci in chromosome 22  
Now working on: INIT POS  
Count statistics have been initialized!  
Count statistics have been initialized!  
Count statistics have been initialized!  
Count statistics have been initialized!  
Count statistics have been initialized!  
Now working on: INIT STATS
```

Once the initial setup is over, The user is prompted to specify QC filters:

```
Indicate the filters and corresponding values(e.g. hwe 10e-5). Available filters are HWE, MAF, MPS(Missing Per sample), MPN(Missing per snp), snp(rsid) (MPS not implemented yet):  
>hwe 1e-5
```

This step is relatively quick and you should see the following response from the server: 

```
in chromosome 21, 2833 snps were deleted and 47390 snps remain  
in chromosome 22, 3103 snps were deleted and 46674 snps remain
```

The initial command prompt reappears. At this point, we can exit to further analyze the data, we can apply more filters or proceed to perform PCA. We will proceed.

```
Looks like everyone is here. Type the stage you want to proceed with? Options: HELP, INIT, QC, PCA, EXIT:  
>pca

Specify filters for PCA pre-processing: (OPTIONS: HWE, MAF, MPS, as before as well as LD. e.g. maf 0.1 LD 50:  
>maf 0.1 ld 50
...
```
