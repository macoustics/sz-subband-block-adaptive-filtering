clear all; close all; clc
SubbandFlag = true;
cubeMeasurement = false;
musicFlag = false;
NumberOfRepeats = 10;
NumberOfCrossoverBands = 3;
addpath('Mex/')
addpath('MexReal/')
% profile on
%% Choose subband settings
% OversamplingFactors = [2, 3/2, 4/3, 5/4, 8/7];
% NumberOfChannelsList = [4, 8, 16, 32, 64, 128];
NumberOfChannelsList = [4];
for LIST_IDX = 1:length(NumberOfChannelsList)
OversamplingFactor = 2;
NumberOfChannels = NumberOfChannelsList(LIST_IDX);
DecimationFactor = round(NumberOfChannels/OversamplingFactor);
PrototypeFilterLength = round(4*NumberOfChannels/(OversamplingFactor-1));

SetupEnvironment = 'Simulation';
EvaluationEnvironment = 'CubeMeasurement';

%% Process input signal(s)
if musicFlag
%     [inSigA,fs] = audioread('../Misc/Deadmau5 - Seeya_48k.wav');
%     [inSigA,fs] = audioread('../Misc/Sigur Ros - Var_48k.wav');
%     [inSigA,fs] = audioread('../Misc/14 Just a Lil Bit_48000kHz_mono.wav');
    [inSigA,fs] = audioread('../Misc/Daft Punk - Give Life Back To Music_48000kHz_mono.wav');
%     [inSigA,fs] = audioread('../Misc/Foo Fighters - The Pretender_48000kHz_mono.wav');
%     [inSigA,fs] = audioread('../Misc/Rage Against the Machine-Wake Up_48000kHz_mono.wav');
    inSigA = inSigA*1;
    [inSigB,fs] = audioread('../Misc/Deadmau5 - Seeya_48k.wav');
    inSigB = inSigB*0;  
    sigLen = 25*fs;
else
    fs = 48e3;
    sigLen = 30*fs;
    rng(1000)
    inSigA = randn(sigLen + fs,1)/10;
    inSigB = zeros(sigLen + fs,1);
    sigLen = 25*fs;
end

% Create x-over networ
lowFreq{1} = 300;
lowFreq{2} = 1500;
lowFreq{3} = 4000;
highFreq{1} = 1500;
highFreq{2} = 4000;
highFreq{3} = 10000;
if cubeMeasurement
    load('../TransferFunctions/Xover_eq_sos_filters.mat');
    for i = 1:NumberOfCrossoverBands
        InputSignalA{i} = sosfilt(sosOut{i},inSigA);
        InputSignalB{i} = sosfilt(sosOut{i},inSigB);
    end
else
    [sosx, ~, ~, ~] = linkwitzrileyhpc([300, 1500, 4000, 10000], fs, 8);
    for i = 1:NumberOfCrossoverBands
        InputSignalA{i} = sosfilt(sosx{i+1},inSigA);
        InputSignalB{i} = sosfilt(sosx{i+1},inSigB);
    end
end

%% Adjust which gains should be active
%%
ZoneAIdx = 15:17;
ZoneBIdx = 21:23;

load(['DecomposedRirs\' SetupEnvironment '_K=' int2str(NumberOfChannels) '_D=' int2str(DecimationFactor) '_Lp=' int2str(PrototypeFilterLength)]);

%% Generate the subband signals
for bandIdx = 1:NumberOfCrossoverBands
    maxIterations = 2;
    SourceReferenceIdx = round(size(DecomposedRirs{bandIdx},3)/2);
    numberOfLoudspeakers = size(DecomposedRirs{bandIdx},3);
    FirTaps = 300;
    ModellingDelay = 150;
    numberOfBlocks = floor(sigLen/DecimationFactor);
    
    %% Perform adaptive processing
    zerosSignal = zeros(length(InputSignalA{bandIdx}),1);
    [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon] = prepareInputsForLeakyNlms(InputSignalA{bandIdx}, zerosSignal, IRs{bandIdx}(:,:,:), FirTaps, SourceReferenceIdx, ZoneAIdx, ZoneBIdx, ModellingDelay);
    
    

    %% Matlab ref
%     tmpAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, FirTaps, brightRir, darkRir, targetFilterSpectra, stepRange, maxIterations, regParameter, signalEpsilon);
%         
%     hopSize = fftSize/2;
%     numberOfBlocks = floor(length(InputSignalA{bandIdx})/hopSize);
%     outputSignals = zeros(numberOfBlocks*hopSize,size(brightRir,2));
%     ooLoudspeakerSignal{bandIdx} = outputSignals;
%     tic
%     for i = 1:numberOfBlocks
%         idx = (i-1)*hopSize + (1:hopSize);
%         tmpSamples = tmpAdaptiveBeamformer.processInputBuffer(InputSignalA{bandIdx}(idx));
%         ooLoudspeakerSignal{bandIdx}(idx,:) = tmpSamples;
%     end
%     disp(['Matlab OO took ' num2str(toc,'%.2f') ' s'])
    %% Subband processing
    % [fullMexLoudspeakerSignal, durationInNanoseconds] = testSubbandAdaptiveFilteringSpeed(PrototypeFilter, InputSignalA{bandIdx}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(ones(NumberOfChannels/2,1)));
    returnStr = '[fullMexLoudspeakerSignal, durationInNanoseconds] = ';
    varStr1 = ['testAdaptiveFilteringMultiMicsSpeedReal'];
%     varStr1 = ['testSubbandAdaptiveFilteringParallel_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
    paramStr = ['(InputSignalA{' int2str(bandIdx) '}, brightRir, darkRir, targetFilterSpectra, signalEpsilon, regParameter, stepRange);'];
%     tic

    switch bandIdx
        case 1
            varStr2 = ['_src9'];
%             [fullMexLoudspeakerSignal, durationInNanoseconds] = testSubbandAdaptiveFilteringSpeed_K8_D7_Lp224_src9(PrototypeFilter, InputSignalA{bandIdx}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(GainFactors{bandIdx}));
        case 2
            varStr2 = '_src11';
%             [fullMexLoudspeakerSignal, durationInNanoseconds] = testSubbandAdaptiveFilteringSpeed_K8_D7_Lp224_src11(PrototypeFilter, InputSignalA{bandIdx}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(GainFactors{bandIdx}));
        case 3
            varStr2 = '_src7';
%             [fullMexLoudspeakerSignal, durationInNanoseconds] = testSubbandAdaptiveFilteringSpeed_K8_D7_Lp224_src7(PrototypeFilter, InputSignalA{bandIdx}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(GainFactors{bandIdx}));
    end
    for repeatIdx = 1:NumberOfRepeats
        tic
        [fullMexLoudspeakerSignal, durationInNanoseconds] = eval([varStr1 varStr2 paramStr]);
        disp(['Matlab duration is ' num2str(toc,'%.2f') ' s. Mex computations took ' num2str(durationInNanoseconds*1e-9,'%.2f') ' s'])
        
        duration{bandIdx}(repeatIdx) = durationInNanoseconds*1e-9;
    end
    %% Estimate result
    max(max(abs(fullMexLoudspeakerSignal)))
    evalBrightRir = squeeze(IRs{bandIdx}(:,ZoneAIdx(2),:));
    evalDarkRir = squeeze(IRs{bandIdx}(:,ZoneBIdx(2),:));
    fullMexBrightPressure{bandIdx} = predictPressureResponse(evalBrightRir, fullMexLoudspeakerSignal);
    fullMexDarkPressure{bandIdx} = predictPressureResponse(evalDarkRir, fullMexLoudspeakerSignal);

end

SaveFolderName = 'TimeResults/';
fname = ['processingDurationReal'];
save([SaveFolderName fname], "duration");

% brightPressure = zeros(length(fullMexBrightPressure{1}),1);
% darkPressure = brightPressure;
% for bandIdx = 1:NumberOfCrossoverBands
%     brightPressure = brightPressure + fullMexBrightPressure{bandIdx};
%     darkPressure = darkPressure + fullMexDarkPressure{bandIdx};
% end
% 
% scale = norm(brightPressure)/sqrt(length(brightPressure));
% 
% figure
% subplot(2,1,1)
% plot(brightPressure/scale)
% hold on; grid on
% legend('full mex')
% title('Bright pressure')
% 
% subplot(2,1,2)
% plot(darkPressure/scale)
% hold on; grid on
% legend('full mex')
% title('Dark pressure')
% 
% PlotDirectivityResults(IRs, ooLoudspeakerSignals)

end
