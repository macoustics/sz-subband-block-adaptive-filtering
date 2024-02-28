function PlotDirectivityResults(EvaluationRoomImpulseResponses, LoudspeakerSignals)
    NumberOfAngles = size(EvaluationRoomImpulseResponses{1},2);
    Pressure = PredictPressure(LoudspeakerSignals, EvaluationRoomImpulseResponses, 1:NumberOfAngles);
    m_Nfft = 2^12;
%     Spectra = malog(mean(EstimateSpectra(Pressure,m_Nfft), 2), 1/24);
    Spectra = EstimateSpectra(Pressure,m_Nfft);

    for idx = 1:NumberOfAngles
        Spectra(:,idx) = malog(Spectra(:,idx), 1/24);
    end

    Fs = 48e3;
    m_PlotFrequencies = 0:Fs/m_Nfft:Fs-Fs/m_Nfft;
    phi = -90:5:90;

    LogSpectra = 10*log10(abs(Spectra));
    LogSpectra = LogSpectra - max(max(LogSpectra));
    figure
    surf(m_PlotFrequencies, phi, LogSpectra','EdgeColor','none'); view(2)
    xlim([100 18e3]); 
    caxis([-30 0])
    set(gca,'xscale','log')

end