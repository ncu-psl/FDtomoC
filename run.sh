#!/bin/sh

# please make sure current dir is in FDtomo/build
pwd

spec_file=${1}

# make1d
./make1d <<< $spec_file

# c2f
./c2f <<< $spec_file

# sphfd
# not yet parallel

# copy nesscressery files to data folder
cp tempvel.pvel ../data/small/TTimes00/VP.mod
cp tempvel.svel ../data/small/TTimes00/VS.mod
cp sphfd ../data/small/TTimes00/sphfd

cd ../data/small/TTimes00
sh runsphfd01.sh
cd ../../../build

# sphfdloc
./sphfdloc <<< $spec_file

# sphrayderv
./sphrayderv <<< $spec_file

# runlsqr
./runlsqr <<< $spec_file

# makenewmod
./makenewmod <<< $spec_file
