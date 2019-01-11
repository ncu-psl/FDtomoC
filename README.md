# FDtomo-C
* 這是C語言版本的 repository
* 測試資料與文件在 [FDtomo](https://github.com/ncu-psl/FDtomo)
* 任何問題請使用 issue

## 其他版本
* [FDtomo-Fortran(Original)](https://github.com/ncu-psl/FDtomo)
* [FDtomo-Python](https://github.com/ncu-psl/FDtomo-Python)

## How to download
```
# clone source code
cd <where you want FDtomo-C to live>
git clone https://github.com/ncu-psl/FDtomo-C.git

# clone test data
git clone https://github.com/ncu-psl/FDtomo.git
cd FDtomo-C
cp -r ../FDtomo/data ./
```

## How to build
```
# GNU Compiler Collection
export CC=gcc

# Intel C Compiler
# export CC=icc

# build
mkdir build
cd build
cmake ../src/
make -j # if you get Internal Compiler Error(ICE), try make without -j
```

## How to run
```
pwd
# please make sure current dir is in FDtomo/build

# make1d
./make1d <<< ../data/small.spec

# c2f
./c2f <<< ../data/small.spec

# sphfd
# not yet parallel
cp sphfd ../data/small/TTimes00/sphfd
cd ../data/small/TTimes00
sh runsphfd01
cd ../../../build

# sphfdloc
./sphfdloc <<< ../data/small.spec

# sphrayderv
./sphrayderv <<< ../data/small.spec

# runlsqr
./runlsqr <<< ../data/small.spec

# makenewmod
./makenewmod <<< ../data/small.spec
```
