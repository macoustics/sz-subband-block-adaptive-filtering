clear all; close all; clc
%% Generate the curves displayed in Fig. 11
% This script, generates the four individual data curves which depict the
% frequency-domain acoustic contrast (i.e., the ratio of mean-square pressures
% in the bright and the dark zone).

task = 'plotDirectivity';
numberOfChannels = 8;
oversamplingFactor = 5/4;
setupEnvironment = 'Simulation';
numberOfRepeats = 1;
musicFlag = false;
switchZoneFlag = false;
decimationFactor = round(numberOfChannels/oversamplingFactor);
prototypeFilterLength = round(4*numberOfChannels/(oversamplingFactor-1));


%% Subband - point source simulation
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, 'simulation', numberOfRepeats, 2, switchZoneFlag);
title('Point-source simulation')

%% Fullrate - variable step-size
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, 'anechoic', numberOfRepeats, 2, switchZoneFlag);
title('Measured responses')