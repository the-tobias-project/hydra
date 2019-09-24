#!/usr/bin/env bash

if ! docker network ls | grep -q hydra-network; then
	docker network create hydra-network
fi

if ! docker network ls | grep -q hydra-redis; then
	docker network create hydra-redis
fi


cp -r ../src ../build/

docker-compose -p hydra -f ./docker-compose.yml --project-directory=../build up -d

docker attach $(docker-compose -p hydra -f ./docker-compose.yml ps -q app)
