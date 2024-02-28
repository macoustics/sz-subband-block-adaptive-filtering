clear all; close all; clc
SubbandFlag = true;
cubeMeasurement = false;
musicFlag = false;
NumberOfRepeats = 10;
NumberOfCrossoverBands = 3;
addpath('Mex/')
addpath('MexParallel/')
addpath('AnalysisFilterBankMex/')
parallelFlag = false;
allChannelsFlag = true;
% profile on
%% Choose subband settings
OversamplingFactorList = [3/2];
% NumberOfChannelsList = [4, 8, 16, 32, 64, 128];
NumberOfChannelsList = [128];
% OversamplingFactorList = [2];
for FACTOR_IDX = 1:length(OversamplingFactorList)
    OversamplingFactor = OversamplingFactorList(FACTOR_IDX);
for LIST_IDX = 1:length(NumberOfChannelsList)

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
n = (0:NumberOfChannels-1)';

SubbandCenterFrequencies = zeros(NumberOfChannels,1);
for i = 0:NumberOfChannels-1
    SubbandCenterFrequencies(i+1) = fs/NumberOfChannels*(i+0.5);
end
SubbandLowerCutoffFrequencies = SubbandCenterFrequencies-0.5*fs/NumberOfChannels;
SubbandLowerCutoffFrequencies(SubbandLowerCutoffFrequencies<0) = 0;
SubbandHigherCutoffFrequencies = SubbandCenterFrequencies+0.5*fs/NumberOfChannels;

oct = 1/3;
for band = 1:NumberOfCrossoverBands
    GainFactors{band} = zeros(NumberOfChannels/2,1);
    for channel = 1:NumberOfChannels/2
        if (SubbandHigherCutoffFrequencies(channel) > lowFreq{band}/2^oct) && (SubbandLowerCutoffFrequencies(channel) < highFreq{band}*2^oct)
            GainFactors{band}(channel) = 1;
        end
    end
end
% keyboard

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
    FirTaps = ceil(FirTaps/DecimationFactor);
    ModellingDelay = ceil(ModellingDelay/DecimationFactor);
    tmpAnalysisFilterbank = analysisFilterbankFast(PrototypeFilter, NumberOfChannels, DecimationFactor);
    numberOfBlocks = floor(sigLen/DecimationFactor);
%     analyzedSignals = zeros(NumberOfChannels/2, numberOfBlocks);
%     for i = 1:numberOfBlocks
%         idx = (i-1)*DecimationFactor + (1:DecimationFactor);
%         tmpSamples = tmpAnalysisFilterbank.processInputBuffer(InputSignalA{bandIdx}(idx));
%         analyzedSignals(:,i) = tmpSamples;
%     end
    varStr1 = ['testAnalysisFilterbankSpeed_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
    paramStr = ['(PrototypeFilter, InputSignalA{' int2str(bandIdx) '});'];
    analyzedSignals = eval([varStr1 paramStr]);
%     analyzedSignals = testAnalysisFilterbankSpeed_K8_D7_Lp224(PrototypeFilter, InputSignalA{bandIdx});
    
    subbandSignalLength = size(analyzedSignals,2);
    
    %% Perform adaptive processing
    refLoudspeakerSubbandSignal = zeros(NumberOfChannels/2, subbandSignalLength, numberOfLoudspeakers);
    ooLoudspeakerSubbandSignal = refLoudspeakerSubbandSignal;
    mexLoudspeakerSubbandSignal = refLoudspeakerSubbandSignal;
    for cIdx = 1:NumberOfChannels/2
        zerosSignal = zeros(size(analyzedSignals,2),1);
        [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon] = prepareInputsForLeakyNlms(analyzedSignals(cIdx,:).', zerosSignal, DecomposedRirs{bandIdx}(:,:,:,cIdx), FirTaps, SourceReferenceIdx, ZoneAIdx, ZoneBIdx, ModellingDelay);
        if cIdx == 1 
            brightRirArray = zeros(size(brightRir,1),size(brightRir,2), size(brightRir,3),NumberOfChannels/2);
            darkRirArray = brightRirArray;
            targetFilterSpectraArray = zeros(size(brightRir,1),size(brightRir,2), NumberOfChannels/2);
            stepRangeArray = zeros(2,NumberOfChannels/2);
            regParameterArray = zeros(NumberOfChannels/2,1);
            signalEpsilonArray = zeros(NumberOfChannels/2,1);
        end
        brightRirArray(:,:,:,cIdx) = brightRir;
        darkRirArray(:,:,:,cIdx) = darkRir;
        targetFilterSpectraArray(:,:,cIdx) = targetFilterSpectra;
        stepRangeArray(:,cIdx) = stepRange;
        regParameterArray(cIdx) = regParameter;
        signalEpsilonArray(cIdx) = signalEpsilon;
    end
    
    %% Subband processing
    % [fullMexLoudspeakerSignal, durationInNanoseconds] = testSubbandAdaptiveFilteringSpeed(PrototypeFilter, InputSignalA{bandIdx}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(ones(NumberOfChannels/2,1)));
    returnStr = '[fullMexLoudspeakerSignal, durationInNanoseconds] = ';
    
    if parallelFlag
        varStr1 = ['testSubbandAdaptiveFilteringMicsParallel_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
    else
        varStr1 = ['testSubbandAdaptiveFilteringMultipleMics_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
    end
        
    if allChannelsFlag
        paramStr = ['(PrototypeFilter, InputSignalA{' int2str(bandIdx) '}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(ones(NumberOfChannels/2,1) ));'];
    else
        paramStr = ['(PrototypeFilter, InputSignalA{' int2str(bandIdx) '}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(GainFactors{' int2str(bandIdx) '}));'];
    end
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
    tmpLoudspeakerSignals{bandIdx} = fullMexLoudspeakerSignal;
    %% Estimate result
    evalBrightRir = squeeze(IRs{bandIdx}(:,ZoneAIdx(3),:));
    evalDarkRir = squeeze(IRs{bandIdx}(:,ZoneBIdx(3),:));
    fullMexBrightPressure{bandIdx} = predictPressureResponse(evalBrightRir, fullMexLoudspeakerSignal);
    fullMexDarkPressure{bandIdx} = predictPressureResponse(evalDarkRir, fullMexLoudspeakerSignal);
end

SaveFolderName = 'TimeResults/';
if parallelFlag
    fname = ['processingDurationParallel_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
elseif allChannelsFlag
    fname = ['processingDurationAllChannels_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
else
    fname = ['processingDuration_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
end
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
% PlotDirectivityResults(IRs, tmpLoudspeakerSignals)

end
end
