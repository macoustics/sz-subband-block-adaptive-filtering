clear all; close all; clc
addpath('PEAQ-master/PQevalAudio');
addpath('PEAQ-master/PQevalAudio/CB');
addpath('PEAQ-master/PQevalAudio/Misc');
addpath('PEAQ-master/PQevalAudio/MOV');
addpath('PEAQ-master/PQevalAudio/Patt');
addpath('PEAQ-master/');

SubbandFlag = true;
cubeMeasurement = false;
musicFlag = true;
NumberOfRepeats = 1;
maxIterations = 2;
NumberOfCrossoverBands = 3;
addpath('Mex/')
addpath('MexParallel/')
addpath('AnalysisFilterBankMex/')
parallelFlag = false;
allChannelsFlag = true;
switchZonesFlag = false;
% profile on
%% Choose subband settings
OversamplingFactorList = [2];
NumberOfChannelsList = [8];
% NumberOfChannelsList = [32];
% OversamplingFactorList = [2];
qualityList;

for oversamplingIdx = 1:1
for filterbankIdx = 7:7%size(filterbankSettings,2)

NumberOfChannels = filterbankSettings{oversamplingIdx,filterbankIdx}.numberOfChannels;
DecimationFactor = filterbankSettings{oversamplingIdx,filterbankIdx}.decimationFactor;
PrototypeFilterLength = filterbankSettings{oversamplingIdx,filterbankIdx}.filterLength;

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
    inSigA = inSigA(31*fs:end);
    inSigA = inSigA*1;
%     [inSigB,fs] = audioread('../Misc/Deadmau5 - Seeya_48k.wav');
    sigLen = 30*fs;
    inSigA = inSigA(1:sigLen);
    inSigB = inSigA*0;  
    
else
    fs = 48e3;
    sigLen = 30*fs;
    rng(1000)
    inSigA = randn(sigLen,1)/10;
    inSigB = zeros(sigLen,1);
%     sigLen = 25*fs;
end

% Create x-over network
lowFreq{1} = 100;
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
    [sosx, ~, ~, ~] = linkwitzrileyhpc([100, 1500, 4000, 10000], fs, 8);
    
    for i = 1:NumberOfCrossoverBands
        InputSignalA{i} = sosfilt(sosx{i+1},inSigA);
        InputSignalB{i} = sosfilt(sosx{i+1},inSigB);
%         if i == 1
%             wGain = 10^(-1/20);
%             [Bw,Aw] = biquad(1500,4,1,11,fs);
%             InputSignalA{i} = filter(Bw,Aw,InputSignalA{i}) * wGain;
%             InputSignalB{i} = filter(Bw,Aw,InputSignalB{i}) * wGain;
%         elseif i == 2
%             [Bm1,Am1] = biquad(2700,-5, 2, 9, fs);
%             [Bm2,Am2] = biquad(2200,4, 2, 9, fs);
%             InputSignalA{i} = filter(Bm1,Am1,InputSignalA{i});
%             InputSignalB{i} = filter(Bm1,Am1,InputSignalB{i});
%             InputSignalA{i} = filter(Bm2,Am2,InputSignalA{i});
%             InputSignalB{i} = filter(Bm2,Am2,InputSignalB{i});
%         else
%             tGain = 10^(6/20);
%             [Bt,At] = biquad(2500,-6, 1/sqrt(2), 10, fs);
%             InputSignalA{i} = filter(Bt,At,InputSignalA{i}) * tGain;
%             InputSignalB{i} = filter(Bt,At,InputSignalB{i}) * tGain;
%         end
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
% ZoneAIdx = 16;
% ZoneBIdx = 22;

load(['DecomposedRirs\Quality_' SetupEnvironment '_K=' int2str(NumberOfChannels) '_D=' int2str(DecimationFactor) '_Lp=' int2str(PrototypeFilterLength)]);

%% Generate the subband signals
for bandIdx = 1:NumberOfCrossoverBands
    SourceReferenceIdx = round(size(DecomposedRirs{bandIdx},3)/2);
    numberOfLoudspeakers = size(DecomposedRirs{bandIdx},3);
    FirTaps = 300;
    ModellingDelay = 150;
    if SubbandFlag
        FirTaps = ceil(FirTaps/DecimationFactor);
        ModellingDelay = ceil(ModellingDelay/DecimationFactor);
        
%         varStr1 = ['testAnalysisFilterbankSpeed_K' int2str(NumberOfChannels) '_D' int2str(DecimationFactor) '_Lp' int2str(PrototypeFilterLength)];
%         paramStr = ['(PrototypeFilter, InputSignalA{' int2str(bandIdx) '});'];
%         analyzedSignals = eval([varStr1 paramStr]);

        tmpAnalysisFilterbank = analysisFilterbankFast(PrototypeFilter, NumberOfChannels, DecimationFactor);
        numberOfBlocks = floor(sigLen/DecimationFactor);
        analyzedSignals = zeros(NumberOfChannels/2, numberOfBlocks);
        for i = 1:numberOfBlocks
            idx = (i-1)*DecimationFactor + (1:DecimationFactor);
            tmpSamples = tmpAnalysisFilterbank.processInputBuffer(InputSignalA{bandIdx}(idx));
            analyzedSignals(:,i) = tmpSamples;
        end
    else
        DecimationFactor = 1;
        NumberOfChannels = 2;
        analyzedSignals = InputSignalA{bandIdx};
    end
    subbandSignalLength = size(analyzedSignals,2);
    
    %% Perform adaptive processing
    refLoudspeakerSubbandSignal = zeros(NumberOfChannels/2, subbandSignalLength, numberOfLoudspeakers);
    
    ooLoudspeakerSubbandSignal = refLoudspeakerSubbandSignal;
    mexLoudspeakerSubbandSignal = refLoudspeakerSubbandSignal;
    firstChannelProcessed = false;
    firstTimeThrough = true;
    for cIdx = 1:NumberOfChannels/2
        if SubbandFlag
            zerosSignal = zeros(size(analyzedSignals,2),1);
            if (GainFactors{bandIdx}(cIdx))
                firstChannelProcessed = true;
                [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon, Wiener1, Wiener2] = prepareInputsForLeakyNlms(analyzedSignals(cIdx,:).', zerosSignal, DecomposedRirs{bandIdx}(:,:,:,cIdx), FirTaps, SourceReferenceIdx, ZoneAIdx, ZoneBIdx, ModellingDelay);

            end
        else
            zerosSignal = zeros(length(analyzedSignals),1);
            firstChannelProcessed = true;
            [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon, Wiener1, Wiener2] = prepareInputsForLeakyNlms(analyzedSignals(:), zerosSignal, IRs{bandIdx}(:,:,:), FirTaps, SourceReferenceIdx, ZoneAIdx, ZoneBIdx, ModellingDelay);
        end
        if firstChannelProcessed && firstTimeThrough
            brightRirArray = zeros(size(brightRir,1),size(brightRir,2), size(brightRir,3),NumberOfChannels/2);
            darkRirArray = brightRirArray;
            targetFilterSpectraArray = zeros(size(brightRir,1),size(brightRir,2), NumberOfChannels/2);
            stepRangeArray = zeros(2,NumberOfChannels/2);
            regParameterArray = zeros(NumberOfChannels/2,1);
            signalEpsilonArray = zeros(NumberOfChannels/2,1);
            hopSize = fftSize/2;
            numberOfBlocks = floor(length(zerosSignal)/hopSize);
            Misalignment{bandIdx} = zeros(numberOfBlocks, NumberOfChannels/2);
            WienerNorm{bandIdx} = zeros(numberOfBlocks, NumberOfChannels/2);
            FilterNorm{bandIdx} = zeros(numberOfBlocks, NumberOfChannels/2);
            StepSizeSave{bandIdx} = zeros(numberOfBlocks, NumberOfChannels/2);
            firstTimeThrough = false;
        end
        brightRirArray(:,:,:,cIdx) = brightRir;
        darkRirArray(:,:,:,cIdx) = darkRir;
        targetFilterSpectraArray(:,:,cIdx) = targetFilterSpectra;
        stepRangeArray(:,cIdx) = stepRange;
        regParameterArray(cIdx) = regParameter;
        signalEpsilonArray(cIdx) = signalEpsilon;

        % Adaptive processing
        tmpAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, FirTaps, brightRir, darkRir, targetFilterSpectra, stepRange, maxIterations, regParameter, signalEpsilon);
        
        
        outputSignals = zeros(numberOfBlocks*hopSize,size(brightRir,2));
        tmpWiener = Wiener1;
        if SubbandFlag
            analysisSignal = analyzedSignals(cIdx,:).';
        else
            analysisSignal = analyzedSignals;
        end
        if (GainFactors{bandIdx}(cIdx)) || (~SubbandFlag)
            tic
            for i = 1:numberOfBlocks
                if i == floor(numberOfBlocks/2)
                    if switchZonesFlag
                        tmpAdaptiveBeamformer.switchZones();
                        tmpWiener = Wiener2;
                    end
                end
                idx = (i-1)*hopSize + (1:hopSize);
                tmpSamples = tmpAdaptiveBeamformer.processInputBuffer(analysisSignal(idx));
                ooLoudspeakerSubbandSignal(cIdx,idx,:) = tmpSamples;
                Misalignment{bandIdx}(i, cIdx) = norm(tmpWiener-tmpAdaptiveBeamformer.getFilterSpectra(),'fro')/norm(tmpWiener,'fro');
                WienerNorm{bandIdx}(i, cIdx) = norm(tmpWiener, 'fro');
                FilterNorm{bandIdx}(i, cIdx) = norm(tmpAdaptiveBeamformer.getFilterSpectra(),'fro');
                StepSizeSave{bandIdx}(i, cIdx) = tmpAdaptiveBeamformer.m_stepSize;
            end
            disp(['Matlab OO took ' num2str(toc,'%.2f') ' s'])
            FinalFilters{bandIdx,cIdx} = tmpAdaptiveBeamformer.getFilterSpectra();
            WienerFilters{bandIdx,cIdx} = tmpWiener;
            PowerSpectrum{bandIdx,cIdx} = tmpAdaptiveBeamformer.m_powerSpectrum;
        end
    end
    
    if SubbandFlag
        numberOfBlocks = size(refLoudspeakerSubbandSignal,2);
        ooLoudspeakerSignal = zeros(numberOfBlocks*DecimationFactor,numberOfLoudspeakers);
        for lIdx = 1:numberOfLoudspeakers
            ooLoudspeakerSignal(:,lIdx) = synthesizeSignal(ooLoudspeakerSubbandSignal(:,:,lIdx), PrototypeFilter, NumberOfChannels, DecimationFactor);
        end
        tmpLoudspeakerSignals{bandIdx} = ooLoudspeakerSignal;
    else
        tmpLoudspeakerSignals{bandIdx} = squeeze(ooLoudspeakerSubbandSignal(1,:,:));
    end
end

% load("DecomposedRirs\CubeMeasurement_K=8_D=4_Lp=32");

ReferenceSignal = zeros(length(InputSignalA{1}),1);
for idx = 1:length(InputSignalA)
    ReferenceSignal = ReferenceSignal + InputSignalA{idx};
end
Pressure = PredictPressure(tmpLoudspeakerSignals, IRs, round(mean(ZoneAIdx)));
[~, irDelay] = max(abs(IRs{3}(:,round(mean(ZoneAIdx)),4)));
filterbankDelay = PrototypeFilterLength-1;

if isnan(filterbankSettings{oversamplingIdx,filterbankIdx}.latency)
    Latency = 1;
else
    Latency = filterbankSettings{oversamplingIdx,filterbankIdx}.latency;
end

ReferenceSignal = ReferenceSignal/sqrt(mean(ReferenceSignal.^2));
Pressure = Pressure(Latency:end-fs)/sqrt(mean(Pressure(Latency:end-fs).^2));
OutputSignalLength = length(Pressure);
ReferenceSignal = ReferenceSignal(1:OutputSignalLength);
Scale = max(abs(ReferenceSignal)) * 1.1;
ReferenceSignal = ReferenceSignal/Scale;
Pressure = Pressure / Scale;

figure
plot(ReferenceSignal);
hold on; grid on
plot(Pressure);
legend('Reference','Pressure')
xlabel('Time [samples]'); ylabel('Magnitude')

audiowrite('ref.wav',ReferenceSignal(3*fs:end),fs);
audiowrite('test.wav',Pressure(3*fs:end),fs);

[PeaqScore, MOV] = PQevalAudio_fn('ref.wav', 'test.wav');
filterbankSettings{oversamplingIdx, filterbankIdx}.peaqScore = PeaqScore;
title(['K = ' int2str(NumberOfChannels) ' D = ' int2str(DecimationFactor) ' Lp = ' int2str(PrototypeFilterLength) ' Peaq = ' num2str(PeaqScore,'%.2f')]);

end
end

% figure
% for fIdx = 1:length(filterbankSettings)
%     plot(filterbankSettings{fIdx}.signalToAliasRatio,filterbankSettings{fIdx}.peaqScore,'xk')
%     hold on; grid on;
% end
% xlabel('SAR [dB]'); ylabel('Peaq ODG')
% title('New calc')
