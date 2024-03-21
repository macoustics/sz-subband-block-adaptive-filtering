clear all; close all; clc
%% Generate the curves displayed in Fig. 10
% This script, generates the four individual data curves which depict the
% frequency-domain acoustic contrast (i.e., the ratio of mean-square pressures
% in the bright and the dark zone).

task = 'plotPressureSpectra';
numberOfChannels = 8;
oversamplingFactor = 5/4;
setupEnvironment = 'Simulation';
evaluationEnvironment = 'anechoic';
numberOfRepeats = 1;
musicFlag = false;
switchZoneFlag = false;
decimationFactor = round(numberOfChannels/oversamplingFactor);
prototypeFilterLength = round(4*numberOfChannels/(oversamplingFactor-1));


%% Subband - variable step-size
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, switchZoneFlag);
title('Subband - variable step-size')


%% Fullrate - variable step-size
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, switchZoneFlag);
title('Fullrate - variable step-size')
