function synthesizedSignal = synthesizeSignal(subbandSignals, prototypeFilter, numberOfSubbandChannels, decimationFactor)

tmpSynthesisFilterbank = synthesisFilterbankFast(prototypeFilter,numberOfSubbandChannels,decimationFactor);
numberOfBlocks = size(subbandSignals,2);
synthesizedSignal = zeros(numberOfBlocks*decimationFactor,1);
for i = 1:numberOfBlocks
    idx = (i-1)*decimationFactor + (1:decimationFactor);
    tmpSamples = tmpSynthesisFilterbank.processInputBuffer(subbandSignals(:,i));
    synthesizedSignal(idx) = tmpSamples;
end