#!/usr/bin/env bash

#change the version number for each new build
docker build -t cmap/base-merino:latest -t cmap/base-merino:v0.0.1 --rm=true .