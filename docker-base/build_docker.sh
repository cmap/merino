#!/usr/bin/env bash

#change the version number for each new build
docker build --platform linux/amd64 -t cmap/base-merino:latest -t cmap/base-merino:v0.0.2 --rm=true .
