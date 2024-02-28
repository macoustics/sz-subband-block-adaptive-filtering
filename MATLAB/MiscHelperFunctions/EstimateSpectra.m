function [AutoSpectra] = EstimateSpectra(Sequences, WindowLength)
% Function for estimating the spectrum of input sequences
%
% =====================================================
% INPUTS:
% -----------------------------------------------------
% Sequences     2D array    Input sequences in time-domain
%                           size = [NumberOfSamples, NumberOfMics]
% WindowLength  Int         Number of frequency bins in the estimated
%                           spectra.
% =====================================================
% OUTPUTS:
% -----------------------------------------------------
% AutoSpectra   2D array    Output spectra in frequency domain
%                           size = [WindowLength, NumberOfMics]
%
m_NumberOfSamples = size(Sequences,1);
m_NumberOfMics = size(Sequences,2);
m_Hopsize = WindowLength;
m_Window = repmat(hann(WindowLength),1,m_NumberOfMics);
m_NumberOfBlocks = floor((m_NumberOfSamples-WindowLength)/m_Hopsize) + 1;
AutoSpectra = zeros(WindowLength, m_NumberOfMics);

for bIdx = 1:m_NumberOfBlocks
    idx = (bIdx-1)*m_Hopsize + (1:WindowLength);
    tmpSequences = Sequences(idx,:) .* m_Window;
    AutoSpectra = AutoSpectra + abs(fft(tmpSequences,WindowLength,1)).^2;
end
AutoSpectra = AutoSpectra/m_NumberOfBlocks;