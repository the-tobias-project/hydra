#!/bin/bash
set -e

echo "Starting SSH ..."
service ssh start
/usr/sbin/sshd

cd src
python -m server

