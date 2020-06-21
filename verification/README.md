# How to verify the output of FDtomoC
1. Check output for some modules(make1d, c2f, sphfd, sphfdloc, makenewmod).
2. We skip checking output of modules sphrayderv, runlsqr, but check output of makenewmod after passing their output to makenewmod.
3. After running multiple iterations of FDtomo and FDtomoC, we check output of sphfdloc only.

# Debugging tools
* Users can use some tools (like diff) to check binaries output of modules(make1d, c2f, sphfd, makenewmod)
* However, when we were in the process of translating FDtomo to FDtomoC, we find that there was some tiny differences between output of FDtomoC & FDtomo.
* And we create this simple tool for us to debug modules make1d, c2f and sphfd.

## How to use
```bash
gcc model_translate.c -o model_translate
./model_translate [output of FDtomo modules] [output of FDtomoC modules] [size of grid]
# example :
# ./model_translate fmod/tempvel cmod/tempvel 52111
```

