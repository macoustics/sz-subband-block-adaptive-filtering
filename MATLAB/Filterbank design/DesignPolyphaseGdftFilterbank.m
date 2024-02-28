function [PrototypeFilter, AnalysisFilters] = DesignPolyphaseGdftFilterbank(gamma, tau, NumberOfChannels, PrototypeFilterLength, DecimationFactor, PlotResultsFlag)
% Function for designing a (near) perfect reconstruction generalized
% Discrete Fourier Transform (complex modulated) filterbank.
%
% ==================================
% INPUTS:
% ----------------------------------
% gamma scalar      Trade-off parameter between stop-band rejection and
%                   perfect reconstruction in prototype filter design.
%                   Default value should be 0.005
% tau   scalar      Parameter for smoothing the iterative prototype filter
%                   design, between iterations.
%                   Default value should be 0.5
% NumberOfChannels  int     Number of subband channels in the filterbank.
%                           Note for real-valued input signals, only the
%                           first NumberOfChannels/2 subbands are required
%                           for reconstruction, due to symmetry.
% PrototypeFilterLength int Number of samples describing the linear-phase
%                           FIR prototype filter.
%                           Default value should be NumberOfChannels * 4 /
%                           (OversamplingFactor-1).
% DecimationFactor int      Decimation of sampling rate in each subband. Is
%                           equal to floor(NumberOfChannels/OversamplingFactor).
% PlotResultsFlag bool      Plots the results of the Filterbank design
%
% ================================
% OUTPUTS:
% --------------------------------
% PrototypeFilter   1D array    Real-valued linear-phase FIR prototype
%                               filter.
% size = [PrototypeFilterLength, 1].
% AnalysisFilters   2D array    Analysis filters, expressed as complex FIR
%                               filters. For this design synthesis filters
%                               are equal to analysis filters.
% size = [PrototypeFilterLength, NumberOfChannels].


%% Check input parameters
if DecimationFactor > NumberOfChannels
    error('DecimationFactor must be lower than the number of filters covering the entire frequency range');
end
if mod(PrototypeFilterLength,NumberOfChannels)
    error('The current implementation assumes the length of the Prototype filter to be integer multiple of the channel number. This is to conveniently create the polyphase representation of the prototype filter.')
end

%% Design prototype filter
maxIterations =  100;
tolerance = 1e-6;
NumberOfDctFrequencies = PrototypeFilterLength; % Perhaps this should be lcm(2*NumberOfChannels,DecimationFactor)

% Initialize first iteration of the prototype filter as a linear phase,
% lowpass filter.
StopbandEdge = 1/DecimationFactor;
PassbandEdge = 1/(NumberOfChannels);
PrototypeFilter0 = firpm(PrototypeFilterLength-1, [0, 0.9*PassbandEdge, PassbandEdge, StopbandEdge, 1.1*StopbandEdge, 1], [1,1,1,1,0,0],[10,2,1]);
PrototypeFilter0 = PrototypeFilter0(:);
SymmetricPrototypeFilter = PrototypeFilter0(1:PrototypeFilterLength/2);

% Design DCT matrix for evaluating stopband energy
StopbandFrequencies = linspace(StopbandEdge,1,NumberOfDctFrequencies)';
k = (0:PrototypeFilterLength-1);
df = StopbandFrequencies(2)-StopbandFrequencies(1);
DctMatrix = cos(pi*StopbandFrequencies*k)/sqrt(pi*df);

% Modulation matrix
ModulationMatrix = zeros(PrototypeFilterLength,NumberOfChannels);
timeIndicees = (0:PrototypeFilterLength-1)';
timeOffset = -(PrototypeFilterLength-1)/2;
for ChannelNumber = 0:NumberOfChannels-1
%     ModulationMatrix(:,ChannelNumber+1) = exp(1j*2*pi*mod((ChannelNumber+0.5)*(timeIndicees + timeOffset)/(NumberOfChannels),1));
    ModulationMatrix(:,ChannelNumber+1) = exp(1j*2*pi*(ChannelNumber+0.5)*(timeIndicees + timeOffset)/(NumberOfChannels));
end

% Symmetry enforsing matrix
% Here we assume PrototypeFilterLength to be even 
S1 = [eye(PrototypeFilterLength/2); fliplr(eye(PrototypeFilterLength/2))];
S2 = zeros(PrototypeFilterLength);
for i = 1:NumberOfChannels
    for j = 1:PrototypeFilterLength/NumberOfChannels
        S2((i-1)*PrototypeFilterLength/NumberOfChannels + j, (j-1)*NumberOfChannels + i) = 1;
    end
end

% Stopband energy
PL = DctMatrix*S1;

% Target response
DesiredPolyphaseResponse = zeros(2*PrototypeFilterLength/NumberOfChannels-1,1);
DesiredPolyphaseResponse(PrototypeFilterLength/NumberOfChannels) = 1/NumberOfChannels;
DesiredStopbandEnergyAtDctFrequencies = zeros(size(DctMatrix,1),1);
TargetResponse = [DesiredPolyphaseResponse; DesiredStopbandEnergyAtDctFrequencies];

% Iterative solution
itr = 1;
PrevSymmetricPrototypeFilter = SymmetricPrototypeFilter;
epsilon = 1;
NormalizedMeanSquareError = zeros(maxIterations,1);
StopbandEnergy = NormalizedMeanSquareError;
AntiIdentityMatrix = fliplr(eye(PrototypeFilterLength/NumberOfChannels));
while (epsilon > tolerance) && (itr < maxIterations)
    % Toeplitz matrix implementing convolution
    V = [];
    PrevPolyphaseFilter = S2*S1*PrevSymmetricPrototypeFilter;
    PPoly = reshape(PrevPolyphaseFilter,PrototypeFilterLength/NumberOfChannels,NumberOfChannels);
    for i = 1:NumberOfChannels
        tmp = toeplitz([PPoly(:,i); zeros(PrototypeFilterLength/NumberOfChannels-1,1)], [PPoly(1,i), zeros(1,PrototypeFilterLength/NumberOfChannels-1)]);
        V = [V, tmp*AntiIdentityMatrix];
    end
    
    A = [V*S2*S1; gamma*PL];
    SymmetricPrototypeFilter = A\TargetResponse;
    SymmetricPrototypeFilter = tau*SymmetricPrototypeFilter + (1-tau)*PrevSymmetricPrototypeFilter;
    NormalizedMeanSquareError(itr) = norm(V*S2*S1*SymmetricPrototypeFilter-DesiredPolyphaseResponse)/norm(DesiredPolyphaseResponse);
    StopbandEnergy(itr) = norm(PL*SymmetricPrototypeFilter);
    epsilon = norm(SymmetricPrototypeFilter-PrevSymmetricPrototypeFilter)/norm(SymmetricPrototypeFilter);
    PrevSymmetricPrototypeFilter = SymmetricPrototypeFilter;
    itr = itr + 1;
end

if PlotResultsFlag
    % Plot cost-function at each iteration
    figure
    plot(20*log10(NormalizedMeanSquareError))
    hold on; grid on
    plot(20*log10(StopbandEnergy))
    legend('Normalized mean square error','Stopband energy')
    xlabel('Iteration number')
    ylabel('[dB]')
    title('Evaluation of the numerically estimated prototype filter')
end


%% Normalize filters
PrototypeFilter = S1*SymmetricPrototypeFilter;
TestAnalysisFilters = zeros(PrototypeFilterLength,NumberOfChannels);
for i = 1:NumberOfChannels
    TestAnalysisFilters(:,i) = ModulationMatrix(:,i).*PrototypeFilter;
end

% Paraunitaryness
NumberOfBlocks = ceil(PrototypeFilterLength/DecimationFactor);
ParaunitaryMatrix = zeros(DecimationFactor,DecimationFactor);
for i = 1:NumberOfBlocks
    if i==NumberOfBlocks
        idx2 = (i-1)*DecimationFactor+1:PrototypeFilterLength;
        Ai = [TestAnalysisFilters(idx2,:); zeros(NumberOfBlocks*DecimationFactor-PrototypeFilterLength,NumberOfChannels)].';
    else
        idx2 = (i-1)*DecimationFactor+1:i*DecimationFactor;
        Ai = TestAnalysisFilters(idx2,:).';
    end
    ParaunitaryMatrix = ParaunitaryMatrix + Ai'*Ai;
end

if PlotResultsFlag
    % Plot results
    Nfft = 2^16;
    PlotFrequencies = 0:2/Nfft:2-2/Nfft;
%     PlotFrequencies = 0:48e3/Nfft:48e3-48e3/Nfft;
    figure
    PrototypeFilterSpectrum = fft(PrototypeFilter,Nfft);
    plot(PlotFrequencies, 20*log10(abs(PrototypeFilterSpectrum)));
    hold on
    for i = 1:NumberOfChannels/2
        TmpSpectrum = fft(TestAnalysisFilters(:,i),Nfft);
        plot(PlotFrequencies, 20*log10(abs(TmpSpectrum)));
        hold on
    end
    grid on
    xlim([0 1])
    xlabel('Normalized Frequency')
    ylabel('Magnitude [dB]')
    title('Filter bank responses')
end

% Normalize the prototype filter to ensure the filter bank will be paraunitary
PrototypeFilter = PrototypeFilter/sqrt(ParaunitaryMatrix(1,1));
NormalizationCoeff = sqrt(ParaunitaryMatrix(1,1));
    

%% Create the filter bank
AnalysisFilters = zeros(PrototypeFilterLength, NumberOfChannels);
for i = 1:NumberOfChannels
    AnalysisFilters(:,i) = ModulationMatrix(:,i).*PrototypeFilter;
end
% Given the design of the filter bank the synthesis filters will be
% identical to the analysis filters
% SynthesisFilters = AnalysisFilters;
end

