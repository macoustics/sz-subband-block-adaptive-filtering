function buildMex(buildTask, decimationFactor, numberOfSubbands, prototypeFilterLength, fftLength, hopSize, numberOfBuffers, filterLength, numberOfLoudspeakers, numberOfMicrophones)

%% Build analysis filterbanks

macroString = sprintf(' -DMY_DECIMATION_FACTOR=%d -DMY_NUMBER_OF_SUBBANDS=%d -DMY_PROTOTYPE_FILTER_LENGTH=%d -DMY_FFT_LENGTH=%d -DMY_HOP_SIZE=%d -DMY_NUMBER_OF_BUFFERS=%d -DMY_FILTER_LENGTH=%d -DMY_NUMBER_OF_LOUDSPEAKERS=%d -DMY_NUMBER_OF_MICROPHONES=%d', decimationFactor, numberOfSubbands, prototypeFilterLength, fftLength, hopSize, numberOfBuffers, filterLength, numberOfLoudspeakers, numberOfMicrophones);
addMexOpts = ' -f mexopts.xml';
linkBoost = ' -IC:\boost_1_82_0\';

switch buildTask
    case 'analysisFilterbank'
        cppFilename = 'testAnalysisFilterbankSpeed';
        outputFolder = '../MATLAB/AdaptiveFilteringTests/AnalysisFilterbankMex/';        
    case 'real'
        cppFilename = 'testAdaptiveFilteringMultiMicsSpeedReal';
        outputFolder = '../MATLAB/AdaptiveFilteringTests/MexReal/';
    case 'mics'
        cppFilename = 'testSubbandAdaptiveFilteringMultipleMics';
        outputFolder = '../MATLAB/AdaptiveFilteringTests/Mex/';
    case 'parallel'
        cppFilename = 'testSubbandAdaptiveFilteringMicsParallel';
        outputFolder = '../MATLAB/AdaptiveFilteringTests/MexParallel/';
end

if strcmp(buildTask,'real')
    outputFilename = cppFilename;
else
    outputFilename = [cppFilename '_K' int2str(numberOfSubbands) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
end
if ~strcmp(buildTask, 'analysisFilterbank')
    outputFilename = [outputFilename '_src' int2str(numberOfLoudspeakers)];
end

eval(['mex' ' -output ' outputFolder outputFilename addMexOpts macroString ' ' cppFilename '.cpp' linkBoost])
end
