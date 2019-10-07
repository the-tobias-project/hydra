# Example GWAS 

In this example, we use HYDRA to run a GWAS on 1K Genome data split between 3 silos using a fake phenotype. The 2504 individual dataset is split into 3 silos N=(900, 900, 704) per silo. For illustration purposes, we have restricted the dataset to 150k snps from chromosomes 20-22.

## Data

Download [`testData.zip`](https://console.cloud.google.com/storage/browser/hydra-example-data) and unzip the folders in the `example` directory. Download the compiled files (ending with .so) and place them in the `src/lib` folder.

If you'd like, you can generate the test dataset using `download1kG.sh` after setting up the container. This is discourages since it involves changing WGS vcf files to plink format and may take a long time. Nevertheless, you can do this using the following command

```bash 
docker exec -it hydra_app_1 "download1kG.sh"
```
If you have downloaded the data, you can ignore that step. You can also use your own data by simply splitting the plink files and placing them in the corresponding data folders.

## Container setup. 

The data has been splitted between 3 centers. Our goal is to simulate three different computers 

We require 2 shell terminals for each center and one shell terminal for the server. Since, we are trying to simulate a distributed setting on a single machine, the setup is slightly different. In particular, we will set up 4 different containers. 1 container for the central hub (hereafter hub) and 1 container for each of the 3 centers. These containers will communicate over a [docker network](https://docs.docker.com/v17.09/engine/userguide/networking/#an-overlay-network-without-swarm-mode) and will only have access to the data that would have been available in the decentralized setting. We will need 2 shell terminals for each center and 1 terminal for the hub.

In the shell screen to be used for server build the container images: 

```bash
bash-3.2$ bash  up.sh
```

This script creates the containers and the network. It also initializes 3 containers for the three centers and attaches the server container to the current shell terminal. You should see the prompt change to `root@hydra`. Alternatively, you can pull the image and initialize the server as follows:

```bash
docker network create hydra-network
docker run --name hydra -p 9001:9001 --hostname hydra -it --network=hydra-network --network-alias=hydra_network --rm apoursh/hydra:0.1 bash 
```

After the container has booted up (either using `up.sh` or pulling from the docker hub), start the server...

```bash 
root@hydra:/app# cd /app/src
root@hydra:/app/src# python -m server --dev 1
 * Serving Flask app "__main__" (lazy loading)
 * Environment: development
 * Debug mode: off
   [INFO ] 2019-09-24 04:52:39,895 /usr/local/lib/python3.6/site-packages/werkzeug/_internal.py                     :: 122 =>  * Running on http://0.0.0.0:9001/ (Press CTRL+C to quit)
```
Head over to `http://localhost:9001/api/ui/` in a browser to monitor and control the server. You should see a page similar to the one in Figure 1 below. 


Next, we will setup all the clients. If `up.sh` was used, these containers have been made and started already so simply attach to each center and run the client...

```bash
bash-3.2$ docker attach Center1
```

If you'd like to pull the containers from the hub, you can use the command: 

```bash
docker network create redis-network
docker run --name Center1 --expose 9001 --hostname hydra -it --network=hydra-network --network-alias=hydra_network --rm  apoursh/hydra:0.1 bash
```

Once the container is setup via either method, run the client

```bash
root@hydra-client:/app# cd src/
root@hydra-client:/app/src# python -m client --name=Center1 --plinkfile=/app/data/dset1  --dev 1
```

You should see a confirmation of registration on the client side and the server terminal. You can also check to see all the registered clients in the UI. 

The last step in the setup is starting up the workers. We need to connect to each container and run a separate worker. 

```bash
root@hydra-client: docker exec -it Center1 "bash"
root@hydra-client: /app# cd src/
C_FORCE_ROOT=1  celery -A celery_worker.celery worker -Q Center1 -n Center1 --concurrency=1
```

Change "Center1" as appropriate for each container. 


## GWAS

Now we are ready to run a GWAS. This can be done through the UI or using `curl` commands. In the UI, click on "Jobs->/tasks/{task_name}" (see Figure 1).

![UI]()

Select "INIT" from the drop-down list of task_names and press "Try it out!" (see Figure 2).

![Starting]()

You should see a confirmation of this command on the server side and soon after on the client side. The worker will provide periodic updates as it works through the initialization steps. Once the process is over you should see an automatically generated figure `QC_pre_filter.png` that can be helpful for filteration (figure 3, left). Proceed by running QC in the same way. You should see a similar post-QC image `QC_post_filter.png`.

![]()

Finally you can proceed with PCA and ASSO.
