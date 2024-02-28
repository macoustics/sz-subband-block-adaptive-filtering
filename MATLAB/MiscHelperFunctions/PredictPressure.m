function [Pressure] = PredictPressure(lsIn, IR, MicIdx)
% INPUTS:
% ======================================
% lsIn     matrix       Prefiltered audio signal for loudspeakers
%                       size = [nSamples, nSrcs]
% IR        Array       Impulse responses from loudspeakers to microphone
%                       positions
%                       size = [nRirTaps, nMics, nSrc]


NumberOfCrossoverBands = length(IR);
si = size(IR{1});
nMics = length(MicIdx);
nRirTaps = si(1);
nSamples = size(lsIn{1},1);
nFFT = nSamples + nRirTaps-1;

Pressure = zeros(nFFT,nMics);
% keyboard

for band = 1:NumberOfCrossoverBands
    SignalSpectrum = fft(lsIn{band},nFFT,1);
    % permute from [nRirTaps, nMics, nSrc] to [nRirTaps, nSrc, nMics] for easier processing
    tmpIR = permute(IR{band},[1,3,2]);
%     keyboard
    for micIdx = 1:nMics
        TransferFunction = fft(tmpIR(:,:,MicIdx(micIdx)),nFFT,1);
        Pressure(:,micIdx) = Pressure(:,micIdx) + ifft( sum(TransferFunction.*SignalSpectrum,2) );
    end
end

end

