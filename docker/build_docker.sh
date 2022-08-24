#!/usr/bin/env bash


#change the version number for each new build
# Note: You must push changes to Github for them to appear in the docker container
docker build --platform linux/amd64 --no-cache -t prismcmap/merino:latest -t prismcmap/merino:v0.1.5 --rm=true .
