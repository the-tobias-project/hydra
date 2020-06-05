#!/usr/bin/env bash

docker-compose -p hydra -f ./build/docker-compose.yml down

if [[ "$(docker network inspect Center --format "{{range .Containers}}T{{end}}")" == "" ]]; then
	docker network rm hydra-network
fi

docker rmi apoursh/hydra:spoke
