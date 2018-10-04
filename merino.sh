#!/usr/bin/env bash

key=$1

case $key in
    assemble)
    ./assemble/assemble.sh $@
    ;;
    card)
    ./card $@
    ;;
    weave)
    ./weave $@
    ;;
    build)
    ./build $@
    ;;
esac
