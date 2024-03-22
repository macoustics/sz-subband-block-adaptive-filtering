# sz-subband-block-adaptive-filtering
Code following the publication of associated paper

Note that the code is currently being refactored, which means that many of the scripts used for running the tests are broken at the moment. The core functions implementing the filterbank and adaptive processing will not change, but the scripts used to generate the results will. This is to limit the amount of manual work required by the user to reproduce the results.

# General description of the repo:
The code presents both MATLAB and C++ implementations of subband block adaptive filters, implemented as object oriented programming.
The adaptive filtering class stores its own state and inputs / outputs buffers of hopSize samples (50% of the FFT size) corresponding to the used overlap-save design used.
The conversion to and from subband signals is performed using the analysis and synthesis filterbank classes. The input (output) to the analysis (synthesis) filterbank is an integer multiple N of decimationFactor (numberOfSubbands/2) samples, while the output (input) of the analysis (synthesis) filterbank is an integer multiple N of numberOfSubband/2 (decimationFactor) samples.

The analysis and synthesis filterbanks rely on the design of a prototype filter. Functions for this design are provided in MATLAB/Filterbank design.

The subband adaptive filtering algorithm requires room impulse responses (or free field responses) which are decomposed into subband approximations. Functions for performing this decomposition are provided in MATLAB/Decompose rirs into subbands. A set of free field point source simulated impulse responses and anechoic measurements of a corresponding loudspeaker array are provided in MATLAB/TransferFunctions.

The MATLAB implementation of the algorithms are located in AdaptiveFilteringTests, while the C++ implementations are in the C++ folder.

# Guide for running the code
Before being able to reproduce all the figures of the article, there are a couple of steps which must be followed as preparation:
 - To create any of the Figures 7 - 15, it is necessary to create the filterbanks and decomposed room impulse responses. This is done by running the script MATLAB/createDecomposedRirs.m
 - Before it is possible to run the speed tests, it is necessary to compile the C++ implementations into .mex files. Once a C++ compiler supporting C++17 has been associated with your MATLAB installation, it should be possible to run the script C++/complieMexFiles.m
 - The Figures 8, 9, and 15 require a 48kHz music file, which should be located in MiscHelperFunctions/music.wav. In the published work, the track Give Life Back to Music by Daft Punk was used in the interval 00:31 to 01:01 (the time code is hard-coded into MATLAB/runAdaptiveProcessing.m and it is expected that the music file is at least 1 min and 1 s long).
 - Figure 15 relies on an implementation of PEAQ. In the published work, we used the implementation by Peter Kabal from https://www.mmsp.ece.mcgill.ca/Documents/Software/index.html. Note that this implementation requires a few modifications to translate the depreciated MATLAB function wavread() to the current audioread(). The main file of the PEAQ evaluation code is expected at the location MATLAB/PEAQ-master/PQevalAudio_fn.m

Once these prerequisites are in place, it should be possible to run the code as the individual functions createFig6.m - createFig15.m located in the MATLAB folder.



# Notes about the C++ implementation:
A couple of things to note about the C++ implementation / compilation to .mex:
 - The C++ code is writting with the assumption of C++17 due to the utilization of unique_ptr and shared_ptr.
 - At the time of writing the code, the GCC compiler included in mingw from the official supported list available for MATLAB2022b was too old to support the required C++ features.
 - A newer version of mingw (including the GCC13.1.0 compiler) was installed on the system and set up to operate with MATLAB.
 - The BOOST libraries v1_82_0 (used for the implementation of a circular buffer) are assumed and should be downloaded seperately and linked during the compilation of the .mex files.
 - The compiler flags are provided in the file C++/mexopts.xml which is called as an argument when compiling to override the default mexopts settings in MATLAB.
 - All C++ datatypes are defined at compile time. This means that seperate compilations of the .mex files must be performed for the desired controlFilterLength, numberOfLoudspeakers, numberOfMicrophones, numberOfBuffers, fftLength, hopSize, decimationFactor, numberOfSubbands, and prototypeFilterLength.
 - Notable optimization flags being used are: -O3 -march=native -mfma -mavx2 and -fopenmp (for the parallel computations)
