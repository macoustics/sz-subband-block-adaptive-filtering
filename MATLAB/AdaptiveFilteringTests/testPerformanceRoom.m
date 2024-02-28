clear all; close all; clc
addpath('../Filterbank/Polyphase filter/')
addpath(pathdef)
SubbandFlag = true;
cubeMeasurement = false;
musicFlag = true;
switchZonesFlag = true;
NumberOfRepeats = 1;
NumberOfCrossoverBands = 3;
% profile on
%% Choose subband settings
% NumberOfChannelsList = [8, 16, 32, 64, 128];
NumberOfChannelsList = [8];
for LIST_IDX = 1:length(NumberOfChannelsList)
OversamplingFactor = 5/4;
NumberOfChannels = NumberOfChannelsList(LIST_IDX);
DecimationFactor = floor(NumberOfChannels/OversamplingFactor);
PrototypeFilterLength = round(4*NumberOfChannels/(OversamplingFactor-1));

SetupEnvironment = 'Simulation';
EvaluationEnvironment = 'CubeMeasurement';

%% Process input signal(s)
if musicFlag
%     [inSigA,fs] = audioread('../Misc/Deadmau5 - Seeya_48k.wav');
%     [inSigA,fs] = audioread('../Misc/Sigur Ros - Var_48k.wav');
%     [inSigA,fs] = audioread('../Misc/14 Just a Lil Bit_48000kHz_mono.wav');
    [inSigA,fs] = audioread('Daft Punk - Give Life Back To Music_48000kHz_mono.wav');
%     [inSigA,fs] = audioread('../Misc/Foo Fighters - The Pretender_48000kHz_mono.wav');
%     [inSigA,fs] = audioread('../Misc/Rage Against the Machine-Wake Up_48000kHz_mono.wav');
    inSigA = inSigA*1;
%     [inSigB,fs] = audioread('../Misc/Deadmau5 - Seeya_48k.wav');
    inSigB = inSigA*0;  
    sigLen = 30*fs;
else
    fs = 48e3;
    sigLen = 30*fs;
    rng(1000)
    inSigA = randn(sigLen + fs,1)/10;
    inSigB = zeros(sigLen + fs,1);
    sigLen = 30*fs;
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
    [sosx, ~, ~, ~] = linkwitzrileyhpc([500, 1500, 4000, 10000], fs, 8);
    
    for i = 1:NumberOfCrossoverBands
        InputSignalA{i} = sosfilt(sosx{i+1},inSigA);
        InputSignalB{i} = sosfilt(sosx{i+1},inSigB);
        if i == 1
            wGain = 10^(-1/20);
            [Bw,Aw] = biquad(1500,4,1,11,fs);
            InputSignalA{i} = filter(Bw,Aw,InputSignalA{i}) * wGain;
            InputSignalB{i} = filter(Bw,Aw,InputSignalB{i}) * wGain;
        elseif i == 2
            [Bm1,Am1] = biquad(2700,-5, 2, 9, fs);
            [Bm2,Am2] = biquad(2200,4, 2, 9, fs);
            InputSignalA{i} = filter(Bm1,Am1,InputSignalA{i});
            InputSignalB{i} = filter(Bm1,Am1,InputSignalB{i});
            InputSignalA{i} = filter(Bm2,Am2,InputSignalA{i});
            InputSignalB{i} = filter(Bm2,Am2,InputSignalB{i});
        else
            tGain = 10^(6/20);
            [Bt,At] = biquad(2500,-6, 1/sqrt(2), 10, fs);
            InputSignalA{i} = filter(Bt,At,InputSignalA{i}) * tGain;
            InputSignalB{i} = filter(Bt,At,InputSignalB{i}) * tGain;
        end
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
        if (SubbandCenterFrequencies(channel) > lowFreq{band}/2^oct) && (SubbandLowerCutoffFrequencies(channel) < highFreq{band}*2^oct)
            GainFactors{band}(channel) = 1;
        end
    end
end

%%
ZoneAIdx = 15:17;
ZoneBIdx = 21:23;
% ZoneAIdx = 16;
% ZoneBIdx = 22;

load(['DecomposedRirs\' SetupEnvironment '_K=' int2str(NumberOfChannels) '_D=' int2str(DecimationFactor) '_Lp=' int2str(PrototypeFilterLength)]);

%% Generate the subband signals
% bandIdx = 1;
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
    analyzedSignals = zeros(NumberOfChannels/2, numberOfBlocks);
    for i = 1:numberOfBlocks
        idx = (i-1)*DecimationFactor + (1:DecimationFactor);
        tmpSamples = tmpAnalysisFilterbank.processInputBuffer(InputSignalA{bandIdx}(idx));
        analyzedSignals(:,i) = tmpSamples;
    end
    
    % analyzedSignals = zeros(NumberOfChannels/2, numberOfBlocks);
    % for channel = 1:NumberOfChannels/2
    %     tmpSignal = filter(AnalysisFilters(:,channel),1,InputSignalA{bandIdx});
    %     tmpSignal = downsample(tmpSignal, DecimationFactor);
    %     analyzedSignals(channel,1:length(tmpSignal)) = tmpSignal;
    % end
    
    
    subbandSignalLength = size(analyzedSignals,2);
    
    %% Perform adaptive processing
    refLoudspeakerSubbandSignal = zeros(NumberOfChannels/2, subbandSignalLength, numberOfLoudspeakers);
    ooLoudspeakerSubbandSignal = refLoudspeakerSubbandSignal;
    first = true;
    for cIdx = 1:NumberOfChannels/2
        if GainFactors{bandIdx}(cIdx)
            zerosSignal = zeros(size(analyzedSignals,2),1);
            [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon] = prepareInputsForLeakyNlms(analyzedSignals(cIdx,:).', zerosSignal, DecomposedRirs{bandIdx}(:,:,:,cIdx), FirTaps, SourceReferenceIdx, ZoneAIdx, ZoneBIdx, ModellingDelay);
            if first == true
                first = false;
    %             keyboard
                brightRirArray = zeros(size(brightRir,1),size(brightRir,2), size(brightRir,3), NumberOfChannels/2);
                darkRirArray = brightRirArray;
                targetFilterSpectraArray = brightRirArray;
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
        
            tmpAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, FirTaps, brightRir, darkRir, targetFilterSpectra, stepRange, maxIterations, regParameter, signalEpsilon);
            
            hopSize = fftSize/2;
            numberOfBlocks = floor(size(analyzedSignals,2)/hopSize);
            outputSignals = zeros(numberOfBlocks*hopSize,size(brightRir,2));
            tic
            for i = 1:numberOfBlocks
                if i == floor(numberOfBlocks/2)
                    if switchZonesFlag
                        tmpAdaptiveBeamformer.switchZones();
                    end
                end
                idx = (i-1)*hopSize + (1:hopSize);
                tmpSamples = tmpAdaptiveBeamformer.processInputBuffer(analyzedSignals(cIdx,idx).');
                ooLoudspeakerSubbandSignal(cIdx,idx,:) = tmpSamples;
            end
            disp(['Matlab OO took ' num2str(toc,'%.2f') ' s'])
        end
    end
    
    %% Synthesize loudspeaker input signals
    numberOfBlocks = size(refLoudspeakerSubbandSignal,2);
    ooLoudspeakerSignal = zeros(numberOfBlocks*DecimationFactor,numberOfLoudspeakers);
    for lIdx = 1:numberOfLoudspeakers
        ooLoudspeakerSignal(:,lIdx) = synthesizeSignal(ooLoudspeakerSubbandSignal(:,:,lIdx), PrototypeFilter, NumberOfChannels, DecimationFactor);
    end
    tmpLoudspeakerSignals{bandIdx} = ooLoudspeakerSignal;
% refLoudspeakerSignal = synthesizeSignalsRef(AnalysisFilters, DecimationFactor, refLoudspeakerSubbandSignal);
% ooLoudspeakerSignal = synthesizeSignalsRef(AnalysisFilters, DecimationFactor, ooLoudspeakerSubbandSignal);
end




%% Estimate Stevie Wonder result
load('DecomposedRirs/NewTruncatedStevieWonderMeasurements600ms.mat');
tmpIRs = cat(2,EvaluationRirZoneA,EvaluationRirZoneB);
MidIdx = [21:29];
SubTwIdx = [9:19];
TwIdx = [1:7];
IRs{1} = tmpIRs(:,:,MidIdx);
IRs{2} = tmpIRs(:,:,SubTwIdx);
IRs{3} = tmpIRs(:,:,TwIdx);

% ZoneAIdx = [2, 5, 8, 11];
% ZoneBIdx = [2, 5, 8, 11]+12;
ZoneAIdx = 1:12;
ZoneBIdx = 13:24;
figure
NumberOfMicrophones = length(ZoneAIdx);
yLimits = [-1 1];
for micIdx = 1:NumberOfMicrophones
    for bandIdx = 1:NumberOfCrossoverBands
        evalBrightRir = squeeze(IRs{bandIdx}(:,ZoneAIdx(micIdx),:));
        evalDarkRir = squeeze(IRs{bandIdx}(:,ZoneBIdx(micIdx),:));
        if bandIdx == 1
            ooBrightPressure = predictPressureResponse(evalBrightRir, tmpLoudspeakerSignals{bandIdx});
            ooDarkPressure = predictPressureResponse(evalDarkRir, tmpLoudspeakerSignals{bandIdx});
        else
            ooBrightPressure = ooBrightPressure + predictPressureResponse(evalBrightRir, tmpLoudspeakerSignals{bandIdx});
            ooDarkPressure = ooDarkPressure + predictPressureResponse(evalDarkRir, tmpLoudspeakerSignals{bandIdx});
        end
    end
    
    
    subplot(NumberOfMicrophones,2,2*(micIdx-1)+1)
    plot(ooBrightPressure(1:end))
    grid on; hold on
    legend('oo')
    title('Bright pressure')
    ylim(yLimits);
    
    subplot(NumberOfMicrophones,2,2*micIdx)
    plot(ooDarkPressure(1:end))
    grid on; hold on
    legend('oo')
    title('Dark pressure')
    ylim(yLimits);
end
sgtitle('Stevie Wonder')

ZoneAIdx = 1:12;
ZoneBIdx = 13:24;
PlotEvaluationResults(IRs, tmpLoudspeakerSignals, ZoneAIdx, ZoneBIdx, InputSignalA);

end
