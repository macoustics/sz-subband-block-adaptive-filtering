function FlopCount = FastSynthesisComplexity(DecimationFactor, NumberOfChannels, PrototypeFilterLength)
% Evaluate the complexity for analyzing DecimationFactor new samples (+
% some old ones) into K/2-subbands.

% Apply conj(D_1) = diag(conj(d_1)) \in \mathbb C^{NumberOfChannels/2 x NumberOfChannels/2}
CostD1 = NumberOfChannels/2*cMult();

% length-NumberOfChannels FFT
CostFft = FourierTransform(NumberOfChannels, false);
% CostFft = 4*NumberOfChannels*log2(NumberOfChannels);

% Add the conjugate NumberOfChannels/2 subband channels (due to real
% input and output of the filterbank)
CostDiscard = 0;

% Apply conj(D_2) = diag(conj(d_1)) \in \mathbb C^{NumberOfChannels x
% NumberOfChannels}
CostD2 = NumberOfChannels*cMult();

% Repeat the samples to the length of the prototype filter
CostRepeat = 0;

% Take DecimationFactor new samples and filter with polyphase prototype
% filter
PolyphaseCost = FastSynthesisPolyphase(DecimationFactor, PrototypeFilterLength);

% Total complexity of synthesizing DecimationFactor new samples from K/2
% subband samples
FlopCount = PolyphaseCost + CostRepeat + CostD2 + CostFft + CostDiscard + CostD1;
end

%% Helper function
function FlopCount = FastSynthesisPolyphase(DecimationFactor, PrototypeFilterLength)
% Multiplication between repeated samples and the prototype filter
CostX1 = PrototypeFilterLength*rMult();

% Add X1 to the output Tapped delay line and remember that the top
% DecimationFactor samples are zero.
CostX2 = (PrototypeFilterLength-DecimationFactor)*rAdd();

% Shift the tapped delay line by DecimationFactor samples to output the
% desired values
CostShift = 0;
FlopCount = CostX1 + CostX2 + CostShift;
end