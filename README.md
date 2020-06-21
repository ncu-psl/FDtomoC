# FDtomoC
* This is the revised and optimized version of [FDtomo] in C language
* Documents and experiment data are collected in [data]
* Any programming issue or logic issue should be tracked in [issues]

[FDtomo]: https://github.com/ncu-psl/FDtomo
[data]: https://github.com/ncu-psl/FDtomo/tree/master/data
[issues]: https://github.com/ncu-psl/FDtomoC/issues

## Variants
* [FDtomoPy] (in Python, under development)

[FDtomoPy]: https://github.com/ncu-psl/FDtomoPy

## How to download
```bash
# clone source code
cd <where you want FDtomoC to live>
git clone https://github.com/ncu-psl/FDtomoC.git

# clone test data
git clone https://github.com/ncu-psl/FDtomo.git
cd FDtomoC
cp -r ../FDtomo/data .
unzip data/cwb_12-15_S/runs_files/arrivals/cwb_12-15_S.zip -d data/cwb_12-15_S/runs_files/arrivals/
```

## How to build
```bash
# GNU Compiler Collection
export CC=gcc

# Intel C Compiler
# export CC=icc

# build
# if you get Internal Compiler Error(ICE), try make without -j
bash build.sh
```

## How to run
```sh
cd build
bash run.sh "../data/small/FDtomo.spec"
```

## Development tree
├─FDtomo   // Fortran version
├─FDtomoC  // Reorganized C version 
├─FDtomo-C-old  //C-Benchmark version
└─FDtomoPy // Python version

