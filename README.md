# HYDRA
Decentralized GWAS

## Setup

The recommended setup requires the following:

* Docker
* Bash
* An internet connection
* The two ports `[9001, 9200]` open on your computer
* The Docker network namespaces `[hydra-network, hydra_redis]` available

### Download the image and run a container
This may be the fastest way, given a good network connection between your machine and
the docker hub:

1. `docker pull edolfostanford/hydra`
2. `docker run --name hydra --expose 9001 --hostname hydra -it edolfostanford/hydra bash`

After running this, you should find yourself inside the docker container with a prompt that looks
like this:

```bash
root@hydra:/app#
```

### Build the image yourself

To run the setup, first navigate to this directory, then run the following:

`bash up.sh`

(or, if the executable bit is active, `./up.sh`)

This will first build your image to include all necessary libraries and runtime requirements, then run the data prep
script located in `testData/`.  This script may take anywhere between 10 minutes to several hours to complete,
depending on network connections and computational resources available.

Once the `up.sh` script has completed, you should find yourself inside the docker container with a prompt that looks
like this:

```bash
root@hydra:/app#
```

Additionally, if you run `docker ps`from the host machine, you should see two new containers currently running - one
for HYDRA itself (e.g. `hydra_app_1`), and one for the Redis instance associated with HYDRA (e.g. `hydra_redis_1`).  
We use Celery to manage client jobs, and Celery uses Redis as its' communication backbone.

To run the experiments on a sample data set, please follow the directions under [data prep](#data-prep-details)

## Running the server, client(s), and worker(s)
You will need at least three terminal sessions of some sort for a minimal run on a single machine.  There are many
strategies for accomplishing this within the same docker session (e.g. tmux, screen, ...), but the recommended method
is to open three separate terminals on the host machine, then connect to the running container with each terminal.  This
allows one to see the state of the entire system with a minimal set of commands.

To connect to a running docker container with a new terminal, assuming the container name is `hydra_app_1`, run:

`docker exec -it hydra_app_1 bash`

#### Running the server
Prerequisites: None

Running the server is the simplest of all:
```bash
cd /app/src
python -m server
```
It should look something like this:
```bash
root@hydra:/app# cd src
root@hydra:/app/src# python -m server
 * Serving Flask app "__main__" (lazy loading)
 * Environment: production
   WARNING: Do not use the development server in a production environment.
   Use a production WSGI server instead.
 * Debug mode: off
[INFO ] 2019-05-08 23:45:17,569 /usr/local/lib/python3.6/site-packages/werkzeug/_internal.py                     :: 122 =>  * Running on http://0.0.0.0:9001/ (Press CTRL+C to quit)
```

From here you can explore the API on a web browser on the same machine by visiting `http://localhost:9001/api/ui/`

The server invoked without any arguments assumes the default configuration inside [`src/lib/settings.py`](src/lib/settings.py).
You can override the defaults by either modifying that file directly, or adding arguments to the command line invocation.
Details can be accessed by calling `python -m server --help`.

#### Running the worker
Prerequisites: None

The worker needs to be associated with the client it's serving - this is done by giving both the same
name.  Here we assume the name is `BioME`, which is within the list of clients inside 
[`src/lib/settings.py::ClientHTTP`](src/lib/settings.py).  Any modifications should be reflected within that class.

```bash
cd /app/src
C_FORCE_ROOT=1 celery -A worker worker -Q BioME -n BioME --concurrency=1
```

Which results in something like this:

```bash
root@hydra:/app/src# C_FORCE_ROOT=1 celery -A worker worker -Q BioME -n BioME --concurrency=1
/usr/local/lib/python3.6/site-packages/celery/platforms.py:796: RuntimeWarning: You're running the worker with superuser privileges: this is
absolutely not recommended!

Please specify a different user using the --uid option.

User information: uid=0 euid=0 gid=0 egid=0

  uid=uid, euid=euid, gid=gid, egid=egid,

 -------------- celery@BioME v4.2.1 (windowlicker)
---- **** -----
--- * ***  * -- Linux-4.9.125-linuxkit-x86_64-with-debian-9.9 2019-05-08 23:45:30
-- * - **** ---
- ** ---------- [config]
- ** ---------- .> app:         cws_queue:0x7faa3183cef0
- ** ---------- .> transport:   redis://hydra_redis:6379//
- ** ---------- .> results:     redis://hydra_redis:6379/
- *** --- * --- .> concurrency: 1 (prefork)
-- ******* ---- .> task events: OFF (enable -E to monitor tasks in this worker)
--- ***** -----
 -------------- [queues]
                .> BioME            exchange=BioME(direct) key=BioME
```
Because we are running inside a sandboxed Docker environment, we can safely ignore the superuser warning.

#### Running the client
Prerequisites: A running server

On startup, the client registers itself with the server - hence the dependency on a running server.  The
server uses this registration to later send out further tasks to each individual client.  On
client SIGTERM or SIGINT, (e.g. `control + c`) it will attempt to unregister itself from the server.

Starting the client also requires some configuration - again we are assuming the client name is `BioME`:

```bash
cd /app/src
python -m client --name=BioME --plinkfile=/app/testData/dset1 
```

This barebones call assumes the defaults inside [`src/lib/settings.py`](src/lib/settings.py), which can be overwritten
either in the file itself, or on the command-line invocation.  Call `python -m client --help` for details.

You may notice that there are three sets of `dset` files produced from the `testData/download1kG.sh`
script - here we are somewhat arbitrarily choosing the first.  If you run multiple clients on the same 
machine, you will need to maintain namespace uniqueness with respect to the `--plinkfile` argument - 
this is because we use this to name an hdf5 file for the client.  


#### Creating a distributed system
An example client docker container is created from the `build/docker-compose.yml` file with service name `client`.  This
client will communicate to the server over the docker network called `hydra_network`, thus allowing one to create two
completely separate containers on the same physical machine.

Of course, outside of a testing environment like the one just described, you can simply specify the IP addresses of
each host.

To connect to the client, first run `docker ps`, then use that output to attach to the client container.  In the
following example, we will assume the server has already been started:

```bash
> docker ps
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                              NAMES
68996b2109cf        hydra_app           "bash"                   7 minutes ago       Up 7 minutes        0.0.0.0:9001->9001/tcp             hydra_app_1
53fd6c270431        hydra_client        "bash"                   7 minutes ago       Up 3 seconds                                           hydra_client_1

> docker attach hydra_client_1
root@hydra-client:/app# curl hydra_network:9001
{
  "detail": "The requested URL was not found on the server. If you entered the URL manually please check your spelling and try again.",
  "status": 404,
  "title": "Not Found",
  "type": "about:blank"
}
```

The `404` message above is success - it shows we can connect to the server!  For testing purposes, you can add
additional services under the `build/docker-compose.yml` file, or create separate docker containers from the base image;
the implementation is left as an exercise to the reader.

## Performing a GWAS
We will start with an overview of the state machine:

```
+----------+          +------------+           +-----------+
|          |          |            |           |           |
|  Server  |          |  Clients   |           |  Server   |
|          +--------->+            +---------->+           |
|  start   |          |  register  |           |  init     |
|          |          |            |           |           |
+----------+          +------------+           +-----+-----+
                                                     |
                                                     |
     +<----------------------------------------------+
     |
     v
+----+-----+          +-----------+
|          |          |           |
|          |          |           |
|    QC    +--------->+    PCA    |
|          |          |           |
|          |          |           |
+----------+          +-----------+
```

Assuming you want `N` clients, once you see that the server has `N` clients registered,
you will need to start the three tasks `[init, qc, pca]`:

```bash
curl http://localhost:9001/api/tasks/INIT -X POST
curl http://localhost:9001/api/tasks/QC -X POST
curl http://localhost:9001/api/tasks/PCA -X POST
```
1\) Initialization: `curl http://localhost:9001/api/tasks/INIT -X POST`

You will need to monitor the server logs to ensure you don't tell the server to start a task before the
preceeding task has completed.  For example, once the initialization task has completed, you should see
a message like the following:

```bash
[INFO ] [...] :: 41 => storing counts
[INFO ] [...] :: 60 => Done getting init reports from clients
[INFO ] [...] :: 61 => Telling clients to store stats
[INFO ] [...] :: 122 => 127.0.0.1 - - [08/May/2019 23:43:18] "POST /api/tasks/INIT/COUNT?client_name=BioME HTTP/1.1" 200 -
```

After which, for each registered client and each chromosome, you should see client responses indicating
they have finished storing their stats, like so:

```bash
[INFO ] [...] :: 45 => [BioME]: Finished with init stats for chrom 20
[INFO ] [...] :: 122 => 172.25.0.1 - - [15/May/2019 22:30:28] "POST /api/clients/BioME/report?status=Finished%20with%20init%20stats%20for%20chrom%2020 HTTP/1.1" 200 -
[INFO ] [...] :: 45 => [BioME]: Finished with init stats for chrom 21
[INFO ] [...] :: 122 => 172.25.0.1 - - [15/May/2019 22:30:29] "POST /api/clients/BioME/report?status=Finished%20with%20init%20stats%20for%20chrom%2021 HTTP/1.1" 200 -
[INFO ] [...] :: 45 => [BioME]: Finished with init stats for chrom 22
[INFO ] [...] :: 122 => 172.25.0.1 - - [15/May/2019 22:30:30] "POST /api/clients/BioME/report?status=Finished%20with%20init%20stats%20for%20chrom%2022 HTTP/1.1" 200 -
```

2\) QC: `curl http://localhost:9001/api/tasks/QC -X POST`

For the QC task, with the `dset1` dataset, you should see the following in the logs of the worker:

```bash
After filtering 20, 67886 snps remain
After filtering 21, 41111 snps remain
After filtering 22, 40898 snps remain

[2019-05-09 01:22:40,433: WARNING/ForkPoolWorker-1] Pefroming QC
[2019-05-09 01:22:40,448: WARNING/ForkPoolWorker-1] After filtering 20, 67886 snps remain
[2019-05-09 01:22:40,509: WARNING/ForkPoolWorker-1] After filtering 21, 41111 snps remain
[2019-05-09 01:22:40,539: WARNING/ForkPoolWorker-1] After filtering 22, 40898 snps remain
[2019-05-09 01:22:40,568: WARNING/ForkPoolWorker-1] Finished reporting counts
```

And in the logs of the server:

```bash
[INFO ] [...] :: 57 => Got task QC/FIN
[INFO ] [...] :: 120 => Done with filtering in QC stage
[INFO ] [...] :: 70 => We can move on
[INFO ] [...] :: 122 => 127.0.0.1 - - [09/May/2019 01:22:40] "POST /api/tasks/QC/FIN?client_name=BioME HTTP/1.1" 200 -
```

3\) PCA: `curl http://localhost:9001/api/tasks/PCA -X POST`

After the PCA task has completed, you should see the following in the logs of the worker:

```bash
[2019-05-15 22:38:31,563: WARNING/ForkPoolWorker-1] Done with LD pruning
[2019-05-15 22:38:36,123: WARNING/ForkPoolWorker-1] Reporting cov: 20_20: (4277, 900) x (900, 4277)
[2019-05-15 22:39:26,014: WARNING/ForkPoolWorker-1] Reporting cov: 21_20: (2806, 900) x (900, 4277)
[2019-05-15 22:39:58,969: WARNING/ForkPoolWorker-1] Reporting cov: 21_21: (2806, 900) x (900, 2806)
[2019-05-15 22:40:25,431: WARNING/ForkPoolWorker-1] Reporting cov: 22_20: (2738, 900) x (900, 4277)
[2019-05-15 22:40:59,882: WARNING/ForkPoolWorker-1] Reporting cov: 22_21: (2738, 900) x (900, 2806)
[2019-05-15 22:41:31,605: WARNING/ForkPoolWorker-1] Reporting cov: 22_22: (2738, 900) x (900, 2738)
[2019-05-15 22:41:54,846: WARNING/ForkPoolWorker-1] Final size will be 9821
```

And in the logs of the server:

```bash
[INFO ] [...] :: 57 => Got task PCA/COV
[INFO ] [...] :: 164 => dealing with 21_21
[INFO ] [...] :: 177 => 21_21
[INFO ] [...] :: 184 => 21_21
[INFO ] [...] :: 186 => Finished storing covariances
[INFO ] [...] :: 122 => 127.0.0.1 - - [09/May/2019 01:33:17] "POST /api/tasks/PCA/COV?client_name=BioME HTTP/1.1" 200 -
```

The remaining snp counts for each chromosome should be exact matches.

As the clients perform their tasks, they may send further subtasks to the server - these 
will show up in the server logs, along with the name of the client that sent the task.
Additionally, once the task is completed, the client will send a status update to the
server, which will again show up in the server logs.

## Data prep details

We will use chromosome 20, 21, and 22 from 1000 Genomes phase 3 for demonstration purposes (n=2504). You
can find the relevant dataset in `testData/`, or download a copy using `testData/download1kG.sh`.  For
the purposes of this demo, the data has been thinned to contain 100k snps, and the individuals are
split into 3 separate datasets.

Running the `testData/download1kG.sh` script will take approximately 20 minutes to one hour to complete, depending
on your network connection and your computing power.  You will also need approximately 45 GB of free disk space.

## Troubleshooting

1.  Issues with `lib.corr`

    You can attempt a recompile from within the `build/` directory with the following:

    `python compiler.py build_ext --inplace`

    Then move the compiled file into `src/lib`.

2.  Issues with opening or creating a file

    This can happen if you attempt to run initialization a second time without first clearing
    the scratch directory.  Try shutting down the server and all clients, clearing the scratch
    directory, then starting them up again.
