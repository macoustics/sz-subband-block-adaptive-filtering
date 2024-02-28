function PlotEvaluationResults(EvaluationRoomImpulseResponses, LoudspeakerSignals, ZoneAIdx, ZoneBIdx, InputSignalA)

%     ReferenceSignal = zeros(length(InputSignalA{1}),1);
%     for idx = 1:length(InputSignalA)
%         ReferenceSignal = ReferenceSignal + InputSignalA{idx};
%     end
    PressureAtoA = PredictPressure(LoudspeakerSignals, EvaluationRoomImpulseResponses, ZoneAIdx);
%     figure
%     plot(ReferenceSignal/sqrt(mean(ReferenceSignal.^2)));
%     hold on; grid on
%     plot(PressureAtoA(:,1)/sqrt(mean(PressureAtoA(:,1).^2)));
%     legend('Reference','Pressure')
%     xlabel('Time [samples]'); ylabel('Magnitude')


    m_Nfft = 2^12;
    SpectraAtoA = EstimateSpectra(PressureAtoA,m_Nfft);
    PressureAtoB = PredictPressure(LoudspeakerSignals, EvaluationRoomImpulseResponses, ZoneBIdx);

    RmsLevel = norm(PressureAtoA,'fro')/sqrt(size(PressureAtoA,1)*size(PressureAtoB,2));
    SpectraAtoB = EstimateSpectra(PressureAtoB, m_Nfft);
    PressureAtoA = PressureAtoA/RmsLevel;
    PressureAtoB = PressureAtoB/RmsLevel;
    
    BlockSize = 1024;
    HopSize = BlockSize/2;
    SignalLength = size(PressureAtoA,1);
    NumberOfBlocks = floor(SignalLength/HopSize) - 1;
    Time = (0:(NumberOfBlocks-1)) * HopSize/48e3;
    TimeContrast = zeros(NumberOfBlocks,1);
    for idx = 1:NumberOfBlocks
        idxStart = (idx-1)*HopSize+1;
        idxEnd = idxStart+BlockSize-1;
        pA = sqrt(mean(mean(PressureAtoA(idxStart:idxEnd,:).^2,2)));
        pB = sqrt(mean(mean(PressureAtoB(idxStart:idxEnd,:).^2,2)));
%         keyboard
        TimeContrast(idx) = 20*log10(pA/pB);
    end

    figure
    plot(Time, TimeContrast)
    xlabel('Time [s]'); ylabel('Contrast [dB]')

    for idx = 1:length(ZoneAIdx)
        SpectraAtoA(:,idx) = malog(SpectraAtoA(:,idx), 1/24);
        SpectraAtoB(:,idx) = malog(SpectraAtoB(:,idx), 1/24);
    end

    Fs = 48e3;
    m_PlotFrequencies = 0:Fs/m_Nfft:Fs-Fs/m_Nfft;

    % Normalization
    LowerFrequency = 500;
    UpperFrequency = 8000;
    [idxL] = find(LowerFrequency<m_PlotFrequencies,1);
    [idxU] = find(UpperFrequency<m_PlotFrequencies,1);
    ReferenceLevel = mean(SpectraAtoA(idxL:idxU));
    figure
    semilogx(m_PlotFrequencies, 10*log10(SpectraAtoA/ReferenceLevel));
    hold on; grid on
    semilogx(m_PlotFrequencies, 10*log10(SpectraAtoB/ReferenceLevel));
    semilogx(m_PlotFrequencies, 10*log10(mean(SpectraAtoA,2)/ReferenceLevel),'b','linewidth',2);
    semilogx(m_PlotFrequencies, 10*log10(mean(SpectraAtoB,2)/ReferenceLevel),'r','linewidth',2);
end