[![DOI](https://zenodo.org/badge/140095242.svg)](https://zenodo.org/badge/latestdoi/140095242)
# Capacity functions

This project contains some simple C functions to evaluate basic information-theoretical quantities in a communication system. 
These functions were designed for optical communications, but they may be applied to any communication system.

All those functions are included in the [capacity_functions.c](capacity_functions.c) file, that is ment to be included in other programs
with the aid of the header file [capacity_functions.h](capacity_functions.h).

Since most of people use [MATLAB](https://www.mathworks.com/), this code includes MEX-functions to call the C functions using MATLAB.

## System requirements and compilation instructions
The code has been designed for 64-bit Linux. Windows-specific instructions are available [here](Windows_Compile.md).

To compile it under Linux, it requires the following software/libraries:
- MATLAB R2018a (or newer versions)
- OpenMP
- GCC
- Make

To compile it, just run `make`.

## Included functions
Full documentation of the function is written using [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

The included functions, for the time being, allow:
- Monte-Carlo evaluation of AWGN Mutual Information (MI) for the AWGN channel for real (PAM) and complex (QAM) modulation formats.
- Calculation of bit-wise Log-Likelihood Ratios (LLRs) for the AWGN channel for real (PAM) and complex (QAM) modulation formats.
- Analytical evaluation of Mutual Information (MI) for the AWGN and AWGN-BICM channel using Gauss-Hermite quadrature for real (PAM) and complex (QAM) modulation 
formats.
- Calculation of bit-wise Log-Likelihood Ratios (LLRs) for the AWGN channel with residual phase noise.

## References
- [F. Kayhan and G. Montorsi, "Constellation Design for Memoryless Phase Noise Channels," in IEEE Transactions on Wireless Communications, vol. 13, no. 5, pp. 2874-2883, May 2014](https://doi.org/10.1109/TWC.2014.040714.130731)
- [A. Alvarado, T. Fehenberger, B. Chen and F. M. J. Willems, "Achievable Information Rates for Fiber Optics: Applications and Computations," in Journal of Lightwave Technology, vol. 36, no. 2, pp. 424-439, Jan 15, 2018.](https://doi.org/10.1109/JLT.2017.2786351)
- [D. Pilori, "Advanced Digital Signal Processing Techniques for High-Speed Optical Links," Ph.D. Thesis, Mar 22, 2019.](https://hdl.handle.net/11583/2729814)

## License
This code is released under [MIT License](https://opensource.org/licenses/MIT).
