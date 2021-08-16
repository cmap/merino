#!/usr/bin/env bash

#change the version number for each new build
docker build -t cmap/merino:latest -t cmap/merino:v0.1.0 --rm=true .
