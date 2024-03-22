clear all; close all; clc
%% Generate the curves displayed in Fig. 13
% This script generates the predicted time-domain contrast between two
% regions in a room. The loudspeaker signals are based on the simulated
% free-field impulse responses. The loudspeaker signals are then convolved
% with the RIRs of the physical loudspeaker array in a real room to two
% microphone array (3 x 4 microphones) positions in the room.

task = 'plotTimeContrast';
numberOfChannels = 8;
oversamplingFactor = 5/4;
setupEnvironment = 'Simulation';
evaluationEnvironment = 'room';
numberOfRepeats = 1;
musicFlag = false;
decimationFactor = round(numberOfChannels/oversamplingFactor);
prototypeFilterLength = round(4*numberOfChannels/(oversamplingFactor-1));


%% Subband - variable step-size
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Subband - variable step-size')
grid on
xlim([0 30]); ylim([-12 12])
