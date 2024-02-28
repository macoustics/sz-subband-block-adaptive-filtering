function FlopCount = FourierTransform(SequenceLength, RealSequenceFlag)
% Flop count for calculating the DFT of an sequence of length
% SequenceLength using a radix-2 Decimation in time algorithm.
if RealSequenceFlag
    % Compensate for spectral buffer-size being half of the time-sequence
    % length
    SequenceLength = 2*(SequenceLength-1); 
    FlopCount = (8/3*SequenceLength*log2(SequenceLength) - 16/9*SequenceLength - 2/9*(-1)^log2(SequenceLength) + 2) * rAdd();
    FlopCount = FlopCount + (4/3*SequenceLength*log2(SequenceLength) - 38/9*SequenceLength + 2/9*(-1)^(log2(SequenceLength)) + 6) * rMult();
else
    FlopCount = (8/3*SequenceLength*log2(SequenceLength) - 16/9*SequenceLength - 2/9*(-1)^log2(SequenceLength) + 2) * rAdd();
    FlopCount = FlopCount + (4/3*SequenceLength*log2(SequenceLength) - 38/9*SequenceLength + 2/9*(-1)^(log2(SequenceLength)) + 6) * rMult();
end
end