#!/bin/bash

# please make sure current dir is in FDtomo/build
pwd

if [ $# -eq 0 ]; then
	echo "Please input spec_file, example:"
	echo ""
	echo "bash run.sh ../data/small/FDtomo.spec"
	echo ""
	exit 1
fi

spec_file=${1}

# make1d
./make1d <<< $spec_file

# c2f
./c2f <<< $spec_file

# sphfd
# not yet parallel

# copy nesscressery files to data folder
cp tempvel.pvel VP.mod
cp tempvel.svel VS.mod
./sphfd <<< $spec_file

# sphfdloc
./sphfdloc <<< $spec_file

# sphrayderv
./sphrayderv <<< $spec_file

# runlsqr
./runlsqr <<< $spec_file

# makenewmod
./makenewmod <<< $spec_file
