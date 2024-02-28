clear all; close all; clc
% Test using a user defined prototype filter

%% Parameters for determining the filterbank
gamma = 0.005;
tau = 0.5;
NumberOfChannels = 32;
OversamplingFactor = 2; % [2/1, 3/2, 5/4]
DecimationFactor = round(NumberOfChannels/OversamplingFactor);
PrototypeFilterLength = round(4*NumberOfChannels/(OversamplingFactor-1));

PlotFlag = true;

[PrototypeFilter, AnalysisFilters] = DesignPolyphaseGdftFilterbank(gamma, tau, NumberOfChannels, PrototypeFilterLength, DecimationFactor, PlotFlag);

%% Analyze the behavior of the filterbank
AnalyzeFilterbankProperties(PrototypeFilterLength, DecimationFactor, PrototypeFilter, AnalysisFilters);
