function FlopCount = FastAnalysisComplexity(DecimationFactor, NumberOfChannels, PrototypeFilterLength)
% Evaluate the complexity for analyzing DecimationFactor new samples (+
% some old ones) into K/2-subbands.

% Take DecimationFactor new samples and filter with polyphase prototype
% filter
PolyphaseCost = FastAnalysisPolyphase(NumberOfChannels, PrototypeFilterLength);

% Apply D_2 = diag(d_2) \in \mathbb C^{NumberOfChannels x NumberOfChannels}
CostD2 = NumberOfChannels*cMult();

% length-NumberOfChannels IFFT
CostIfft = FourierTransform(NumberOfChannels, false);

% Discard the redundant NumberOfChannels/2 subband channels (due to real
% input and output of the filterbank)
CostDiscard = 0;

% Apply D_1 = diag(d_1) \in \mathbb C^{NumberOfChannels/2 x
% NumberOfChannels/2}
CostD1 = NumberOfChannels/2*cMult();

% Total complexity of analyzing DecimationFactor new samples into K/2
% subband samples
FlopCount = PolyphaseCost + CostD2 + CostIfft + CostDiscard + CostD1;
end

%% Helper function
function FlopCount = FastAnalysisPolyphase(NumberOfChannels, PrototypeFilterLength)
% Pprime = P.*TDL
% P = [p(0), ..., p(PrototypeFilterLength-1)].';
% TDL: Tapped delay line of input signal
% TDL \in mathbb R^{PrototypeFilterLength \times 1}
CostX1 = PrototypeFilterLength*rMult();

% Sum into the K subbands
ChannelMultiplier = floor(PrototypeFilterLength/NumberOfChannels);
FilterLengthModChannels = mod(PrototypeFilterLength,NumberOfChannels);
CostX2 = ((ChannelMultiplier-1)*NumberOfChannels + FilterLengthModChannels) * rAdd();

FlopCount = CostX1 + CostX2;
end
