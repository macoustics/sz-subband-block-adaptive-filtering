function AnalyzeFilterbankProperties(PrototypeFilterLength, DecimationFactor, PrototypeFilter, AnalysisFilters)
NumberOfChannels = size(AnalysisFilters,2);
%% Evaluate properties of filter bank
%% Orthogonality of overlapping blocks
NumberOfBlocks = ceil(PrototypeFilterLength/DecimationFactor);
OrthogonalityBetweenBlocks = zeros(DecimationFactor, DecimationFactor*(NumberOfBlocks-1));
for shift = 1:NumberOfBlocks-1
    for i = shift+1:NumberOfBlocks
        if i==NumberOfBlocks
            idx2 = (i-1)*DecimationFactor+1:PrototypeFilterLength;
            Ai2 = [AnalysisFilters(idx2,:); zeros(NumberOfBlocks*DecimationFactor-PrototypeFilterLength,NumberOfChannels)].';
        else
            idx2 = (i-1)*DecimationFactor+1:i*DecimationFactor;
            Ai2 = AnalysisFilters(idx2,:).';
        end
        idx1 = (i-shift-1)*DecimationFactor+1:(i-shift)*DecimationFactor;
        Ai1 = AnalysisFilters(idx1,:).';
        OrthogonalityBetweenBlocks(:,(shift-1)*DecimationFactor + (1:DecimationFactor)) = OrthogonalityBetweenBlocks(:,(shift-1)*DecimationFactor + (1:DecimationFactor)) + Ai2'*Ai1;
    end
end
figure
imagesc(20*log10(abs(OrthogonalityBetweenBlocks)))
caxis([-60 0])
title('Evaluation of the orthogonality between blocks (should be all zeros)')


%% Check for para-unitary
%% Orthogonality of overlapping blocks
NumberOfBlocks = ceil(PrototypeFilterLength/DecimationFactor);
ParaunitaryMatrix = zeros(DecimationFactor,DecimationFactor);
for i = 1:NumberOfBlocks
    if i==NumberOfBlocks
        idx2 = (i-1)*DecimationFactor+1:PrototypeFilterLength;
        Ai = [AnalysisFilters(idx2,:); zeros(NumberOfBlocks*DecimationFactor-PrototypeFilterLength,NumberOfChannels)].';
    else
        idx2 = (i-1)*DecimationFactor+1:i*DecimationFactor;
        Ai = AnalysisFilters(idx2,:).';
    end
    ParaunitaryMatrix = ParaunitaryMatrix + Ai'*Ai;
end
figure
imagesc(20*log10(abs(ParaunitaryMatrix)))
caxis([-60 0])
title('Evaluation of paraunitary matrix (should be an identity matrix)')



%% Signal-to-alias ratio
figure
Nfft = 2^16;
PlotFrequencies = 0:2/Nfft:2-2/Nfft;
PrototypeFilterSpectrum = fft(PrototypeFilter,Nfft);
plot(PlotFrequencies, 20*log10(abs(PrototypeFilterSpectrum)));
hold on; grid on
plot([1,1]*1/DecimationFactor, [-200, 5],'--k','linewidth',2) % Stopband edge
plot([1,1]*1/(NumberOfChannels),[-200, 5],'--r','linewidth',2) % Passband edge
idx1 = round(Nfft/(2*DecimationFactor));
SignalEnergy = sum(abs(PrototypeFilterSpectrum(1:idx1)).^2);
StopbandEnergy = sum(abs(PrototypeFilterSpectrum(idx1+1:Nfft/2+1)).^2);
SignalToAliasRatio = 10*log10(SignalEnergy/StopbandEnergy);
title(['Prototype filter response | SAR = ' num2str(SignalToAliasRatio,'%.2f') ' dB'])
xlim([0 0.5])
xlabel('Normalized frequency')
ylabel('Magnitude [dB]')
legend('Prototype filter','1/DecimationFactor','1/NumberOfChannels')

%% Reconstruction error
RefImpulseResponse = zeros(2*PrototypeFilterLength-1,1);
RefImpulseResponse(PrototypeFilterLength) = 1;
SynthesisFilters = AnalysisFilters;
ChannelImpulseResponses = zeros(2*PrototypeFilterLength-1,NumberOfChannels);
for i = 1:NumberOfChannels
    ChannelImpulseResponses(:,i) = conv(AnalysisFilters(:,i),SynthesisFilters(:,i),'full');
end
FilterbankImpulseResponse = real(sum(ChannelImpulseResponses,2))/DecimationFactor;
ReconstructionError = 20*log10(norm(FilterbankImpulseResponse-RefImpulseResponse)/norm(RefImpulseResponse));

figure
plot(real(ChannelImpulseResponses),'r')
hold on
plot(FilterbankImpulseResponse,'k','LineWidth',2)
title(['Channel transfer impulse responses (black line should be an impulse) | Err = ' num2str(ReconstructionError,'%.2f') ' dB'])
xlabel('Time [samples]')
ylabel('Amplitude')





end