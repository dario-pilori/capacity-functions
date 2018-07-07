# Windows compilation instructions
This project has been developed with Linux in mind. However, it is technically possibile to compile it under Windows using any MATLAB-supported C compilers.

This guide uses the [MSYS2](http://www.msys2.org/) software, which is based on the [MinGW-w64](https://mingw-w64.org/doku.php) compiler.

1. Install and update 64-bit [MSYS2](http://www.msys2.org/) following the instructions on the website.
2. Install 64-bit [MinGW-w64](https://mingw-w64.org/doku.php) in MSYS2 `pacman -S mingw-w64-x86_64-gcc`
3. Export the location of the MEX executable: `export PATH="${PATH}:/c/Program Files/MATLAB/R2018a/bin/win64"`
4. Export the location of MinGW-w64 for MEX: `export MW_MINGW64_LOC=/mingw64`
5. Compile this project using `make` (like Linux)