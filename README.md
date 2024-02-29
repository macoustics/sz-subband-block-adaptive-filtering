# sz-subband-block-adaptive-filtering
Code following the publication of associated paper

The code presents both MATLAB and C++ implementations of subband block adaptive filters, implemented as object oriented programming.
The adaptive filtering class stores its own state and inputs / outputs buffers of hopSize samples (50% of the FFT size) corresponding to the used overlap-save design used.
The conversion to and from subband signals is performed using the analysis and synthesis filterbank classes. The input (output) to the analysis (synthesis) filterbank is an integer multiple N of decimationFactor (numberOfSubbands/2) samples, while the output (input) of the analysis (synthesis) filterbank is an integer multiple N of numberOfSubband/2 (decimationFactor) samples.

The analysis and synthesis filterbanks rely on the design of a prototype filter. Functions for this design are provided in MATLAB/Filterbank design.

The subband adaptive filtering algorithm requires room impulse responses (or free field responses) which are decomposed into subband approximations. Functions for performing this decomposition are provided in MATLAB/Decompose rirs into subbands. A set of free field point source simulated impulse responses and anechoic measurements of a corresponding loudspeaker array are provided in MATLAB/TransferFunctions.

The MATLAB implementation of the algorithms are located in AdaptiveFilteringTests, while the C++ implementations are in the C++ folder.

A couple of things to note about the C++ implementation / compilation to .mex:
 - The C++ code is writting with the assumption of C++17 due to the utilization of unique_ptr and shared_ptr.
 - At the time of writing the code, the GCC compiler included in mingw from the official supported list available for MATLAB2022b was too old to support the required C++ features.
 - A newer version of mingw (including the GCC13.1.0 compiler) was installed on the system and set up to operate with MATLAB.
 - The BOOST libraries v1_82_0 (used for the implementation of a circular buffer) are assumed and should be downloaded seperately and linked during the compilation of the .mex files.
 - The compiler flags are provided in the file C++/mexopts.xml which is called as an argument when compiling to override the default mexopts settings in MATLAB.
 - All C++ datatypes are defined at compile time. This means that seperate compilations of the .mex files must be performed for the desired controlFilterLength, numberOfLoudspeakers, numberOfMicrophones, numberOfBuffers, fftLength, hopSize, decimationFactor, numberOfSubbands, and prototypeFilterLength.
 - Notable optimization flags being used are: -O3 -march=native -mfma -mavx2 and -fopenmp (for the parallel computations)
