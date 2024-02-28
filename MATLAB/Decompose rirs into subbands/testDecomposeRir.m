clear all; close all; clc
addpath('../TransferFunctions')
addpath('../Filterbank design')

%% Parameters for determining the filterbank
gamma = 0.005;
tau = 0.5;
RIRs = 'RoomMeasurement'; % [Simulation, CubeMeasurement, RoomMeasurement]
NumberOfChannels = 32;
OversamplingFactor = 2; % [2/1, 3/2, 5/4]
DecimationFactor = round(NumberOfChannels/OversamplingFactor);

PrototypeFilterLength = round(4*NumberOfChannels/(OversamplingFactor-1));
% PrototypeFilterLength = round(4*NumberOfChannels);
PlotFlag = false;

[PrototypeFilter, AnalysisFilters] = DesignPolyphaseGdftFilterbank(gamma, tau, NumberOfChannels, PrototypeFilterLength, DecimationFactor, PlotFlag);

switch lower(RIRs)
    case 'simulation'
        load('TransferFunctions/rawFreeFieldSimulations_cropped.mat');
        IRs{1} = IRs{1}(:,:,3:11);
    case 'cubemeasurement'
        load('TransferFunctions/Prototype2_Setup1_FreeFieldResponses')
        tmpIRs = resp;
        MidIdx = [21:29];
        SubTwIdx = [9:19];
        TwIdx = [1:7];
        IRs{1} = tmpIRs(:,:,MidIdx);
        IRs{2} = tmpIRs(:,:,SubTwIdx);
        IRs{3} = tmpIRs(:,:,TwIdx);
    case 'roommeasurement'
        load('TransferFunctions/NewTruncatedStevieWonderMeasurements600ms');
        tmpIRs = cat(2,SetupRirZoneA,SetupRirZoneB);
        MidIdx = [21:29];
        SubTwIdx = [9:19];
        TwIdx = [1:7];
        IRs{1} = tmpIRs(:,:,MidIdx);
        IRs{2} = tmpIRs(:,:,SubTwIdx);
        IRs{3} = tmpIRs(:,:,TwIdx);
    otherwise
        error('No valid RIRs have been selected');
end

%% Design the filter bank
[DecomposedRirs] = DecomposeRirs(AnalysisFilters, DecimationFactor, IRs);

%% Save the results
Foldername = 'DecomposedRirs/';
Filename = [RIRs '_K=' int2str(NumberOfChannels) '_D=' int2str(DecimationFactor) '_Lp=' int2str(PrototypeFilterLength)];

save([Foldername Filename], "PrototypeFilter", "AnalysisFilters", "DecimationFactor", "DecomposedRirs", "IRs");