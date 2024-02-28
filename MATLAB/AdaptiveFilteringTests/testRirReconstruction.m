clear all; close all; clc

addpath('AnalysisFilterBankMex/')
MicIdx = 19;

% profile on
%% Choose subband settings
OversamplingFactorList = [5/4];
NumberOfChannelsList = [8];
% NumberOfChannelsList = [32];
% OversamplingFactorList = [2];
Fs = 48e3;
signalLength = Fs;
qualityList;

%% LatencyTable
latTable{1,1} = 29;
latTable{1,2} = 45;
latTable{1,3} = 61;
latTable{1,4} = 77;
latTable{1,5} = 93;
latTable{1,6} = 109;
latTable{1,7} = 125;

latTable{2,1} = 28;
latTable{2,2} = 44;
latTable{2,3} = 60;
latTable{2,4} = 76;
latTable{2,5} = 92;
latTable{2,6} = 108;
latTable{2,7} = 124;

latTable{3,1} = 27;
latTable{3,2} = 43;
latTable{3,3} = 59;
latTable{3,4} = 75;
latTable{3,5} = 91;
latTable{3,6} = 107;
latTable{3,7} = 123;


%%
for oversamplingIdx = 1:3
for filterbankIdx = 1:7%size(filterbankSettings,2)

NumberOfChannels = filterbankSettings{oversamplingIdx,filterbankIdx}.numberOfChannels;
DecimationFactor = filterbankSettings{oversamplingIdx,filterbankIdx}.decimationFactor;
PrototypeFilterLength = filterbankSettings{oversamplingIdx,filterbankIdx}.filterLength;

SetupEnvironment = 'Simulation';
load(['DecomposedRirs\Quality_' SetupEnvironment '_K=' int2str(NumberOfChannels) '_D=' int2str(DecimationFactor) '_Lp=' int2str(PrototypeFilterLength)]);

%% Process input signal(s)
figure
for bandIdx = 1:3
    InputSignal = zeros(signalLength,1);
    InputSignal(1) = 1;
    SrcIdx = round(size(IRs{bandIdx},3)/2);
    RefResponse = filter(IRs{bandIdx}(:,MicIdx, SrcIdx),1,InputSignal);
    
    tmpAnalysisFilterbank = analysisFilterbankFast(PrototypeFilter, NumberOfChannels, DecimationFactor);
    numberOfBlocks = floor(signalLength/DecimationFactor);
    analyzedSignals = zeros(NumberOfChannels/2, numberOfBlocks);
    for i = 1:numberOfBlocks
        idx = (i-1)*DecimationFactor + (1:DecimationFactor);
        tmpSamples = tmpAnalysisFilterbank.processInputBuffer(InputSignal(idx));
        analyzedSignals(:,i) = tmpSamples;
    end
    
    tmpAnalyzedSignal = zeros(size(analyzedSignals,1), size(analyzedSignals,2));
    for cIdx = 1:NumberOfChannels/2
        tmpRir = DecomposedRirs{bandIdx}(:,MicIdx,SrcIdx,cIdx);
        tmpAnalyzedSignal(cIdx,:) = filter(tmpRir,1,analyzedSignals(cIdx,:));
    end
    
    testResponse = synthesizeSignal(tmpAnalyzedSignal, PrototypeFilter, NumberOfChannels, DecimationFactor);
    
    testResponse = testResponse(latTable{oversamplingIdx,filterbankIdx}:end);
    responseLength = length(testResponse);
    RefResponse = RefResponse(1:responseLength);
    reconstructionError = 20*log10(norm(testResponse-RefResponse)/norm(RefResponse));
    
    subplot(3,1,bandIdx)
    plot(RefResponse)
    hold on; grid on;
    plot(testResponse)
    title(['Reconstruction error = ' num2str(reconstructionError,'%.2f') ' dB'])
    xlabel('Time [samples]')
    xlim([1 1000])
end


    

end
end
