clear all; close all; clc

addpath('Filterbank design\');
addpath('Decompose rirs into subbands\')

%% Generic filterbank design parameters
gamma = 0.005;
tau = 0.5;
plotFlag = false;

%% Create decomposed RIRs for Fig. 7, 8, 9, 10, 11, 13, and 14
% numberOfSubbandChannelsList = 2.^(2:7);
numberOfSubbandChannelsList = 2.^(3);
oversamplingFactorList = [5/4];

% Load simulated RIRs
load('TransferFunctions/rawFreeFieldSimulations_cropped.mat');
IRs{1} = IRs{1}(:,:,3:11);

% Specify output folder
outputFolder  = 'Decompose rirs into subbands/DecomposedRirs/';

for cIdx = 1:length(numberOfSubbandChannelsList)
    for oIdx = 1:length(oversamplingFactorList)
        numberOfChannels = numberOfSubbandChannelsList(cIdx);
        oversamplingFactor = oversamplingFactorList(oIdx);
        decimationFactor = round(numberOfChannels/oversamplingFactor);
        prototypeFilterLength = round(4*numberOfChannels/(oversamplingFactor-1));
        [prototypeFilter, analysisFilters] = DesignPolyphaseGdftFilterbank(gamma, tau, numberOfChannels, prototypeFilterLength, decimationFactor, plotFlag);
        [decomposedRirs] = DecomposeRirs(analysisFilters, decimationFactor, IRs);
        save([outputFolder 'Simulation_K=' int2str(numberOfChannels) '_D=' int2str(decimationFactor) '_Lp=' int2str(prototypeFilterLength)], 'decimationFactor', 'decomposedRirs', 'IRs', 'prototypeFilter');
    end
end

%% Create decomposed RIRs for Fig. 15
