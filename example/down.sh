#!/usr/bin/env bash

docker-compose -p hydra -f ./docker-compose.yml down

if [[ "$(docker network inspect hydra-client --format "{{range .Containers}}T{{end}}")" == "" ]]; then
	docker network rm hydra-network
fi

docker rmi hydra_app
