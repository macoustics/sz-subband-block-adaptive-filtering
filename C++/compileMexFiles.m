clear all; close all; clc
% Introduce a pause between compilations to ensure minGW has finished
% before starting the next job
pauseInSeconds = 1;

%% %% Generate list of parameters to use for the build of the .mex files
numberOfMicrophones = 3; % (in each zone)
loudspeakerDriverList = [7, 9, 11];
% numberOfChannels = [4,8,16,32,64,128].';
% oversamplingFactors = [2, 3/2, 5/4];
numberOfChannels = [8].';
oversamplingFactors = [5/4];
% loudspeakerDriverList = [7];
decimationFactors = round(numberOfChannels*(1./oversamplingFactors));

filterLength = 300;
rirLength = 700;
fullRateFftLength = 2^nextpow2(2*(filterLength + rirLength -1)-1);
fullRateHopSize = fullRateFftLength/2;


subsystemFilterLengths = ceil(filterLength./decimationFactors);
prototypeFilterLengths = round(4*numberOfChannels*(1./(oversamplingFactors-1)));
subsampledConvolutionLengths = ceil((rirLength + prototypeFilterLengths-1)./decimationFactors); 
subsystemResponseLengths = subsampledConvolutionLengths - ceil(prototypeFilterLengths./decimationFactors) + 1;

convolutionLengths = subsystemResponseLengths + subsystemFilterLengths - 1;

nfft = 2.^nextpow2(2*convolutionLengths-1);
hopSizes = nfft/2; % The output only takes the control filter
blockSizes = decimationFactors.*hopSizes;
samplingRate = 48e3;
signalDuration = 30;
numberOfBlocks = floor(samplingRate*signalDuration./blockSizes);
fullRateNumberOfBlocks = floor(samplingRate*signalDuration/fullRateHopSize);

%% Compile analysis filterbanks
for cIdx = 1:length(numberOfChannels)
    for oIdx = 1:length(oversamplingFactors)
        disp(['Compiling analysis filterbank cIdx=' int2str(cIdx) '/' int2str(length(numberOfChannels)) ' oIdx=' int2str(oIdx) '/' int2str(length(oversamplingFactors))])
        buildMex('analysisFilterbank', decimationFactors(cIdx,oIdx), numberOfChannels(cIdx), prototypeFilterLengths(cIdx,oIdx), nfft(cIdx,oIdx), hopSizes(cIdx,oIdx), numberOfBlocks, subsystemFilterLengths(cIdx,oIdx), 0, numberOfMicrophones );
        pause(pauseInSeconds)
    end
end

%% Compile subband adaptive filtering
for cIdx = 1:length(numberOfChannels)
    for oIdx = 1:length(oversamplingFactors)
        for lIdx = 1:length(loudspeakerDriverList)
            disp(['Compiling subband adaptive filtering cIdx=' int2str(cIdx) '/' int2str(length(numberOfChannels)) ' oIdx=' int2str(oIdx) '/' int2str(length(oversamplingFactors))])
            buildMex('mics', decimationFactors(cIdx,oIdx), numberOfChannels(cIdx), prototypeFilterLengths(cIdx,oIdx), nfft(cIdx,oIdx), hopSizes(cIdx,oIdx), numberOfBlocks, subsystemFilterLengths(cIdx,oIdx), loudspeakerDriverList(lIdx), numberOfMicrophones );
            pause(pauseInSeconds)
            buildMex('parallel', decimationFactors(cIdx,oIdx), numberOfChannels(cIdx), prototypeFilterLengths(cIdx,oIdx), nfft(cIdx,oIdx), hopSizes(cIdx,oIdx), numberOfBlocks, subsystemFilterLengths(cIdx,oIdx), loudspeakerDriverList(lIdx), numberOfMicrophones );
            pause(pauseInSeconds)
        end
    end
end

%% Compile full-rate adaptive filtering
disp(['Compiling full-rate adaptive filtering'])
for lIdx = 1:length(loudspeakerDriverList)
    buildMex('real', 0, 0, 0, fullRateFftLength, fullRateHopSize, fullRateNumberOfBlocks, filterLength, loudspeakerDriverList(lIdx), numberOfMicrophones );
    pause(pauseInSeconds)
end

