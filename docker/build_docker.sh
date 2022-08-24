#!/usr/bin/env bash

#change the version number for each new build
docker build --platform linux/amd64 --no-cache -t prismcmap/merino:latest -t prismcmap/merino:v0.1.5 --rm=true .
