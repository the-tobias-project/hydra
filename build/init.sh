#!/bin/bash
set -e

echo "Starting SSH ..."
service ssh start

cd src
python -m server

