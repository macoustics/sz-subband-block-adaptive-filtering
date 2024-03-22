clear all; close all; clc
%% Generate the curves displayed in Fig. 14
% This script, generates the two data curves which depict the
% frequency-domain acoustic contrast (i.e., the ratio of mean-square pressures
% in the bright and the dark zone). 

task = 'calculatePeaq';
setupEnvironment = 'Simulation';
evaluationEnvironment = 'simulation';
numberOfRepeats = 1;
musicFlag = true;
switchZoneFlag = false;

%% Subband - variable step-size
load('MiscHelperFunctions\filterbankSettingsVsPrototypeFilterLength.mat','filterbankSettings');
for oIdx = 1:size(filterbankSettings,1)
    for lIdx = 1:size(filterbankSettings,2)
        numberOfChannels = filterbankSettings{oIdx,lIdx}.numberOfChannels;
        decimationFactor = filterbankSettings{oIdx,lIdx}.decimationFactor;
        prototypeFilterLength = filterbankSettings{oIdx,lIdx}.filterLength;
        latency = filterbankSettings{oIdx,lIdx}.latency;
        peaqScore = runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, switchZoneFlag, latency);
        filterbankSettings{oIdx,lIdx}.peaqScore = peaqScore;
    end
end


figure
for idx1 = 1:size(filterbankSettings,1)
    tmpX = [];
    tmpY = [];
    for idx2 = 1:size(filterbankSettings,2)
        tmpX = [tmpX; filterbankSettings{idx1,idx2}.filterLength];
        tmpY = [tmpY; filterbankSettings{idx1,idx2}.peaqScore];
    end
    plot(tmpX,tmpY)
    hold on; grid on
end
xlabel('FilterLength'); ylabel('PEAQ ODG')
legend('O = 2','O = 3/2','O = 5/4')
ylim([-4 0.1])
set(gca,'XTick',[32:16:128])
set(gca,'YTick',[-4:0])