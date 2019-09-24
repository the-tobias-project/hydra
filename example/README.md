# Example GWAS 

In this example, we use HYDRA to run a GWAS on 1K Genome data split between 3 silos using a fake phenotype. 

## Data

Download the [data](https://console.cloud.google.com/storage/browser/hydra-example-data) and place it in a folder named testData and unzip the files.

## Container setup. 

We will split the data between 3 centers. We require 2 shell terminals for each center and one shell terminal for the server. In the server shell terminal build the container images: 

```bash
bash-3.2$ bash  up.sh
```

This script creates the containers and initializes 3 containers for the three centers. At the end of the script, it attaches the server container. Start the server...

```bash 
root@hydra:/app# cd /app/src
root@hydra:/app/src# python -m server --dev 1
 * Serving Flask app "__main__" (lazy loading)
 * Environment: development
 * Debug mode: off
   [INFO ] 2019-09-24 04:52:39,895 /usr/local/lib/python3.6/site-packages/werkzeug/_internal.py                     :: 122 =>  * Running on http://0.0.0.0:9001/ (Press CTRL+C to quit)
```
At this point, you can head over to `http://localhost:9001/api/ui/` in a browser to monitor and control the server.


Next, we will setup all the clients. Attach to each server and run the client...

```bash
bash-3.2$ docker attach Center1
root@hydra-client:/app# cd src/
root@hydra-client:/app/src# python -m client --name=Center1 --plinkfile=/app/data/dset1  --dev 1
```

You should see a confirmation of registration on the client side and the server terminal. You can also check to see all the registered clients in the UI. 

The last step in the setup is starting up the workers

