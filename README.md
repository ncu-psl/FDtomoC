# FDtomo-C
* 這是C語言版本的 repository
* 測試資料與文件在 [FDtomo](https://github.com/ncu-psl/FDtomo)
* 任何問題請使用 issue

## 其他版本
* [FDtomo-Fortran(Original)](https://github.com/ncu-psl/FDtomo)
* [FDtomo-Python](https://github.com/ncu-psl/FDtomo-Python)

## How to download
```sh
# clone source code
cd <where you want FDtomo-C to live>
git clone https://github.com/ncu-psl/FDtomo-C.git

# clone test data
git clone https://github.com/ncu-psl/FDtomo.git
cd FDtomo-C
cp -r ../FDtomo/data ./
```

## How to build
```sh
# GNU Compiler Collection
export CC=gcc

# Intel C Compiler
# export CC=icc

# build
# if you get Internal Compiler Error(ICE), try make without -j
sh build.sh
```

## How to run
```sh
cd build
sh run.sh "../data/small/FDtomo.spec"
```
