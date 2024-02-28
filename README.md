# sz-subband-block-adaptive-filtering
Code following the publication of associated paper

A couple of things to note about the C++ implementation / compilation to .mex:
 - The C++ code is writting with the assumption of C++17 due to the utilization of unique_ptr and shared_ptr.
 - At the time of writing the code, the GCC compiler included in mingw from the official supported list available for MATLAB2022b was too old to support the required C++ features.
 - A newer version of mingw (including the GCC13.1.0 compiler) was installed on the system and set up to operate with MATLAB.
 - The BOOST libraries v1_82_0 are assumed and should be downloaded seperately and linked during the compilation of the .mex files.
 - The compiler flags are provided in the file C++/mexopts.xml which is called as an argument when compiling to override the default mexopts settings in MATLAB.
 - All C++ datatypes are defined at compile time. This means that seperate compilations of the .mex files must be performed for the desired controlFilterLength, numberOfLoudspeakers, numberOfMicrophones, numberOfBuffers, fftLength, hopSize, decimationFactor, numberOfSubbands, and prototypeFilterLength.
 - Notable optimization flags being used are: -O3 -march=native -mfma -mavx2 and -fopenmp (for the parallel computations)
