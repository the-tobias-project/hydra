#!/usr/bin/env bash

if ! docker network ls | grep -q hydra-network; then
	docker network create hydra-network
fi

cp testData/download1kG.sh build/

docker-compose -p hydra -f ./build/docker-compose.yml up -d

docker attach $(docker-compose -p hydra -f ./build/docker-compose.yml ps -q app)
