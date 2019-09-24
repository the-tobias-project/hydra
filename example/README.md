# Example GWAS 

In this example, we use HYDRA to run a GWAS on 1K Genome data split between 3 silos using a fake phenotype. The 2504 individual dataset is split into 3 silos N=(900, 900, 704) per silo. For illustration purposes, we have restricted the dataset to 150k snps from chromosomes 20-22.

## Data

Create the `testData` subdirectory (`mkdir testData`). Download the Center specific [data](https://console.cloud.google.com/storage/browser/hydra-example-data) and unzip the folders in `testData`. Download the compiled files (ending with .so) and place them in the src/lib folder.

You can generate the test dataset using `download1kG.sh` after setting up the container using the following command

```bash 
docker exec -it hydra_app_1 "download1kG.sh"
```
Note that processing the whole genome files is slow. You can also use your own data by simply splitting the plink files and placing them in the corresponding data folders.

## Container setup. 

The data has been splitted between 3 centers. We require 2 shell terminals for each center and one shell terminal for the server. In the shell screent to be used for server build the container images: 

```bash
bash-3.2$ bash  up.sh
```

This script creates the containers and initializes 3 containers for the three centers. After starting up the containers the script attaches the server container to the current shell terminal and you should see the prompt change to `root@hydra`. Start the server...

```bash 
root@hydra:/app# cd /app/src
root@hydra:/app/src# python -m server --dev 1
 * Serving Flask app "__main__" (lazy loading)
 * Environment: development
 * Debug mode: off
   [INFO ] 2019-09-24 04:52:39,895 /usr/local/lib/python3.6/site-packages/werkzeug/_internal.py                     :: 122 =>  * Running on http://0.0.0.0:9001/ (Press CTRL+C to quit)
```
At this point, you can head over to `http://localhost:9001/api/ui/` in a browser to monitor and control the server.


Next, we will setup all the clients. On a different shell terminal, attach to each server and run the client...

```bash
bash-3.2$ docker attach Center1
root@hydra-client:/app# cd src/
root@hydra-client:/app/src# python -m client --name=Center1 --plinkfile=/app/data/dset1  --dev 1
```

You should see a confirmation of registration on the client side and the server terminal. You can also check to see all the registered clients in the UI. 

The last step in the setup is starting up the workers. We need to connect to each container and run a separate worker. 

```bash
root@hydradocker exec -it Center1 "bash"
root@hydra-client:/app# cd src/
C_FORCE_ROOT=1  celery -A celery_worker.celery worker -Q Center1 -n Center1 --concurrency=1
```

Change "Center1" as appropriate for each container. 


## GWAS

Now we are ready to run a GWAS. This can be done through the UI or using `curl` commands. In the UI, click on "tasks->/tasks/{task_name}". Select "INIT" from the drop-down list of task_names and press "Try it out!". You should see a confirmation of this command on the server side and soon after on the client side. The worker will provide periodic updates as it works through the initialization steps. Once the process is over follow up with QC, PCA and ASSO in the same way. 



