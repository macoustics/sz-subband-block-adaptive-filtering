clear all; close all; clc
%% Generate the curves displayed in Fig. 9
% This script, generates the four individual data curves which depict the
% time-domain acoustic contrast (i.e., the ratio of mean-square pressures
% in the bright and the dark zone).

task = 'plotTimeContrast';
numberOfChannels = 8;
oversamplingFactor = 5/4;
setupEnvironment = 'Simulation';
evaluationEnvironment = 'anechoic';
numberOfRepeats = 1;
musicFlag = true;
decimationFactor = round(numberOfChannels/oversamplingFactor);
prototypeFilterLength = round(4*numberOfChannels/(oversamplingFactor-1));


%% Subband - fixed step-size
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Subband - fixed step-size')

%% Subband - variable step-size
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Subband - variable step-size')

%% Fullrate - fixed step-size
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Fullrate - fixed step-size')

%% Fullrate - variable step-size
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Fullrate - variable step-size')
