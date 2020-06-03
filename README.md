# HyDRA
Software for decentralized GWAS

## Tutorial  

The easiest way to familiarize yourself with the pipeline is following this short [example](https://github.com/apoursh/HyDRA/tree/master/example) file. The example will walk through the setup and a mock GWAS using publicly available 1K genomes data.


## Setup 

HyDRA is consisted of two separate pieces. The server/hub and the client(s)/spoke(s). In this section, we explain how to setup either section. Most users will be interested in the spoke setup. 

For both sections the recommended setup requires the following: 

* Docker (Version 2.0 or above) 
* Bash
* An internet connection
* The two ports `[9001, 9200]` open on your computer
* The Docker network namespaces `[hydra-network, hydra_redis]` available

### Spoke 
The Docker image can be downloaded or build from scratch. To download the image run:
```
docker pull apoursh/hydra:spoke
```
then run `bash up.sh`. This command will build the container if it has not been downloaded. 

<!--
Server side: 
1. `docker pull apoursh/hydra:spoke`
2. `docker run --name hydra -p 9001:9001 --hostname hydra -it apoursh/hydra:0.1 bash`
-->

You should find yourself inside the docker container with a prompt that looks like this:

```bash
root@hydra:/app#
```
<!--
### Build the image yourself 

Instead of downloading the image, you can build the image from scratch. To run the setup first download the binary files (ending with `.so`) from [here](https://console.cloud.google.com/storage/browser/hydra-example-data) and place them in the `src/lib` directory. Then navigate to the software-home directory and run the following:

`bash up.sh`

(or, if the executable bit is active, `./up.sh`)

This will build your image to include all necessary libraries and runtime requirements. Once the `up.sh` script has completed, you should find yourself inside the docker container with a prompt that looks like this:

```bash
root@hydra:/app#
```
-->

If you run `docker ps`from the host machine, you should see two new containers currently running - one
for HyDRA itself (e.g. `Center`), and one for the Redis instance associated with HyDRA (e.g. `hydra_redis_1`). 
    
We use Celery to manage client jobs, and Celery uses Redis as its communication backbone. Next we will setup the working queue.

We will use the current window for running the spoke but we will need another window to run the redis worker on. Open another terminal and connect to the container using the command:

```
docker exec -it  Center bash
```
change directory to the src directory `cd src` and start the worker: 

```bash
cd /app/src
C_FORCE_ROOT=1 celery -A worker worker -Q Center -n Center --concurrency=1
```

Which results in something like this:

```bash
root@hydra:/app/src# C_FORCE_ROOT=1 celery -A worker worker -Q Center1 -n Center1 --concurrency=1
/usr/local/lib/python3.6/site-packages/celery/platforms.py:796: RuntimeWarning: You're running the worker with superuser privileges: this is
absolutely not recommended!

Please specify a different user using the --uid option.

User information: uid=0 euid=0 gid=0 egid=0

  uid=uid, euid=euid, gid=gid, egid=egid,

 -------------- celery@Center1 v4.2.1 (windowlicker)
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
                .> Center1            exchange=Center(direct) key=Center1
```
Because we are running inside a sandboxed Docker environment, we can safely ignore the superuser warning.

Now we are ready to run the spoke (client)

On startup, the client registers itself with the server. On
client SIGTERM or SIGINT, (e.g. `control + c`) it will attempt to unregister itself from the server.

Starting the client also requires some configuration - For all examples, we assume the client name is `Center1`:

```bash
cd /app/src
python -m client --name=Center1 --plinkfile=/app/testData/dset1 
```

This barebones call assumes the defaults inside [`src/lib/settings.py`](src/lib/settings.py), which can be overwritten
either by directly modifying the file, or on the command-line invocation.  Call `python -m client --help` for details.


### Hub

The study organizer (/organization) will need to setup the hub. Below are the instructions for setting up the hub using Azure App Services. 

First, set the default parameters in [`src/lib/settings.py`](src/lib/settings.py). The `external_host` variable corresponds to the name of the app website. 

Next, build the app image file and push to dockerhub:

```
docker build . --file build/Dockerfile --tag <dockerhub/name:hub>
docker push <dockerhub/name:hub>
```

Our image can be found at `apoursh/hydra:hub`

To implement the image as an app, log into your Azure account and open an Azure Cloud Shell terminal. 

I assume a resource group `RG` has been made for this app. 

```
az appservice plan create --name hydraAppService --resource-group RG --sku <B1> --is-linux
az webapp create --resource-group RG --plan hydraAppService  --name hydraapp --deployment-container-image-name <dockerhub/name:hub>
az webapp config appsettings set --resource-group RG --name hydraapp --settings WEBSITES_PORT=9001
```
The first line creates an app service plan. In this example I have used a B1 plan but you can scale that to your needs. The second line deploys the app and the last line configures the app service to use the port 9001 for the website (this was set in `ServerHTTP` class of [`src/lib/settings.py`](src/lib/settings.py). the names `hydraAppService, hydraapp` can be changed as desired. 

It will take a few minutes for the container to download. When ready the server can be accessed at `https://hydraapp.azurewebsites.net/api/ui` (replace hydraapp with the name you had chosen). 

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
+----+-----+          +-----------+          +-----------+
|          |          |           |          |           |
|          |          |           |          |           |
|    QC    +--------->+    PCA    +--------->+    ASSO   |
|          |          |           |          |           |
|          |          |           |          |           |
+----------+          +-----------+          +-----------+
```



<!---
## Below is depricated
#### Creating a distributed system
An example client docker container is created from the `build/docker-compose.yml` file with service name `client`.  This
client will communicate to the server over the docker network called `hydra_network`, thus allowing one to create two
completely separate containers on the same physical machine.

Of course, outside of a testing environment like the one just described, you can simply specify the IP addresses of
each host.

To connect to the client, first run `docker ps`, then use that output to attach to the client container.  In the
following example, we will assume the server has already been started:
-->
<!--
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

Assuming you want `N` clients, once you see that the server has `N` clients registered,
you will need to sequentially start the tasks `[init, qc, pca, asso]`:

```bash
curl http://localhost:9001/api/tasks/INIT -X POST
curl http://localhost:9001/api/tasks/QC -X POST
curl http://localhost:9001/api/tasks/PCA -X POST
curl http://localhost:9001/api/tasks/ASSO -X POST
```
These tasks can also be started using the UI. 


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


As the clients perform their tasks, they may send further subtasks to the server - these 
will show up in the server logs, along with the name of the client that sent the task.
Additionally, once the task is completed, the client will send a status update to the
server, which will again show up in the server logs.


## Parameter costumization

The default parameters are specified in `src/lib/settings.py` under `Thresholds` class. To change these defaults, the user can directly modify the file (only the server side file is relevant) or specify the desired threshold using the same keywords before executing a task (e.g. `{"QC_HWE": "1e-6", "QC_MAF": ".05"}`). The parameters specified using the UI are not case sensitive and take precedence over the parameters specified in the `settings.py` file. Any extra parameters is ignored. 

Currently the uptions are:

* QC_hwe
* PCA_maf
* PCA_ld_window
* PCA_ld_threshold
* PCA_pcs
* ASSO_pcs



## Troubleshooting

1.  Issues with opening or creating a file

    This can happen if you attempt to run initialization a second time without first clearing
    the scratch directory or if the program exists unexpectedly. If you would like to re-initialize 
    the experiment, try shutting down the server, all clients and workers. Then remove the 
    scratch directory, and restart.

    If the process was unexpectedly stopped and exited, you can clear the file pointer flag using 
    `h5clear -s scratch/filename.h5py`.

## Other resources 

Note that while we only transfer aggregate data, we do this multiple times as a result, just like 
publication of summary statistics for meta-analysis (and even more so), this work is susceptible to
Homer-like attacks by the server. If the server cannot be trusted encryption based methods such as 
[this work](https://github.com/hhcho/secure-gwas) are more appropriate. 


1.  Issues with `lib.corr`

    You can attempt a recompile from within the `build/` directory with the following:

    `python compiler.py build_ext \-\-inplace`

    Then move the compiled file into `src/lib`.
-->
