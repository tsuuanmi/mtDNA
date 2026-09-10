#!/bin/bash
export LC_ALL=C
ROOTDIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
export PATH=$ROOTDIR/mfold-3.6/bin:$ROOTDIR/gs10.02.0/bin:$PATH

"$ROOTDIR/mfold-3.6/bin/mfold" "${@:1}"
