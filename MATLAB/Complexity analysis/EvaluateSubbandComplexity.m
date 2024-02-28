function [CostFullband, CostSubband, CostSubbandActiveChannel] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, NumberOfChannels, DecimationFactor, PrototypeFilterLength, NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequency, fs)
%% Input variables
BlockSize = 2^nextpow2(2*(RirLength+LoudspeakerFilterLength-1)-1);
HopSize = BlockSize/2;

%% Filterbank variables
DecimatedRirLength = ceil((PrototypeFilterLength+RirLength-1)/DecimationFactor) - ceil(PrototypeFilterLength/DecimationFactor) + 1;
DecimatedLoudspeakerFilterLength = ceil(LoudspeakerFilterLength/DecimationFactor);
DecimatedBlockSize = 2^nextpow2(2*(DecimatedRirLength+DecimatedLoudspeakerFilterLength-1)-1);
DecimatedHopSize = DecimatedBlockSize/2;

%% Complexity of fullband filtering per sample
NumberOfCrossoverBands = length(NumberOfLoudspeakers);
CostFullband = 0;
for band = 1:NumberOfCrossoverBands
    tmp = GradientDescentComplexity(BlockSize, NumberOfLoudspeakers(band), NumberOfMicrophones, NumberOfGradientSteps, NumberOfLineSearchSteps, true);
    CostFullband = CostFullband + tmp/HopSize;
end

CostSubbandAdaptiveFilter = 0;
CostSubbandAdaptiveFilterRecursive = 0;
CostSubbandSynthesis = 0;
tmp = DecimatedHopSize*FastAnalysisComplexity(DecimationFactor, NumberOfChannels, PrototypeFilterLength);
CostSubbandAnalysis = tmp/(DecimatedHopSize*DecimationFactor);
for band = 1:NumberOfCrossoverBands
    tmp = DecimatedHopSize*FastSynthesisComplexity(DecimationFactor, NumberOfChannels, PrototypeFilterLength);
    CostSubbandSynthesis = CostSubbandSynthesis + tmp/(DecimatedHopSize*DecimationFactor) * NumberOfLoudspeakers(band);

    tmp = GradientDescentComplexity(DecimatedBlockSize, NumberOfLoudspeakers(band), NumberOfMicrophones, NumberOfGradientSteps, NumberOfLineSearchSteps, false);
    CostSubbandChannelAdaptation = tmp/(DecimatedHopSize*DecimationFactor);
    CostSubbandAdaptiveFilter = CostSubbandAdaptiveFilter + CostSubbandChannelAdaptation*(NumberOfChannels/2);
end
CostSubband(1) = CostSubbandAdaptiveFilter;
CostSubband(2) = CostSubbandAnalysis;
CostSubband(3) = CostSubbandSynthesis;

%% Complexity of reduced bandwidth filtering per sample
CostSubbandActiveChannelAdaptiveFilter = 0;
oct = 1/3;
for band = 1:NumberOfCrossoverBands
    CenterFrequency = zeros(NumberOfChannels,1);
    for i = 0:NumberOfChannels/2-1
        CenterFrequency(i+1) = fs/NumberOfChannels*(i+0.5);
    end
    LowFrequencyLimit = CenterFrequency-0.5*fs/NumberOfChannels;
    LowFrequencyLimit(LowFrequencyLimit<0) = 0;
    HighFrequencyLimit = CenterFrequency+0.5*fs/NumberOfChannels;
    ActiveChannelIndicator = zeros(NumberOfChannels/2,1);
    for channel = 1:NumberOfChannels/2
        if (HighFrequencyLimit(channel) > SourceLowCutoffFrequency(band)/2^oct) && (LowFrequencyLimit(channel) < SourceHighCutoffFrequency(band)*2^oct)
            ActiveChannelIndicator(channel) = 1;
        end
    end
    tmp = GradientDescentComplexity(DecimatedBlockSize, NumberOfLoudspeakers(band), NumberOfMicrophones, NumberOfGradientSteps, NumberOfLineSearchSteps, false);
    CostSubbandChannelAdaptation = tmp/(DecimatedHopSize*DecimationFactor);
    CostSubbandActiveChannelAdaptiveFilter = CostSubbandActiveChannelAdaptiveFilter + CostSubbandChannelAdaptation*sum(ActiveChannelIndicator);
end
CostSubbandActiveChannel(1) = CostSubbandActiveChannelAdaptiveFilter;
CostSubbandActiveChannel(2) = CostSubbandAnalysis;
CostSubbandActiveChannel(3) = CostSubbandSynthesis;

disp(['K=' int2str(NumberOfChannels), ' D=' int2str(DecimationFactor), ' Lp=' int2str(PrototypeFilterLength), ' BlkSize=' int2str(DecimatedBlockSize), ' AF Fullband=' num2str(CostFullband,'%.2e') ' SB=' num2str(sum(CostSubband),'%.2e') ' Synth=' num2str(CostSubbandSynthesis,'%.2e') ' AF=' num2str(CostSubbandActiveChannelAdaptiveFilter,'%.2e')])
