#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

if [ ! -d $DIR/bamtools/src ]; then
    git submodule init
    git submodule update
fi

#mkdir -p $DIR/deps
#if [ ! -d $DIR/deps/bamtools ]; then
#    git clone https://github.com/pezmaster31/bamtools.git $DIR/deps/bamtools
#    cd $DIR/deps/bamtools
#    git reset --hard cbca090607f2154aa19c20ca778c3d067384b606
#    git apply $DIR/patches/bamtools.patch
#    cd $DIR
#fi

if [ ! -d $DIR/deps/bamtools/build ]; then
    mkdir $DIR/deps/bamtools/build
    cd $DIR/deps/bamtools/build
    cmake .. && make
fi



