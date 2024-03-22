function [varargout] = runAdaptiveProcessing(task, algorithm, numberOfChannels, decimationFactor, prototypeFilterLength, musicFlag, setupEnvironment, evaluationEnvironment, numberOfRepeats, maxIterations, switchZonesFlag, varargin)

numberOfCrossoverBands = 3;
addpath('AdaptiveFilteringTests\Mex\')
addpath('AdaptiveFilteringTests\MexParallel\')
addpath('AdaptiveFilteringTests\MexReal\')
addpath('AdaptiveFilteringTests\AnalysisFilterBankMex\')
addpath('AdaptiveFilteringTests\')
addpath('MiscHelperFunctions\')


%% Choose control algorithm
switch algorithm
    case 'matlabFullrate'
        subbandFlag = false;
        mexFlag = false;
    case 'matlabSubband'
        subbandFlag = true;
        allChannelsFlag = false;
        mexFlag = false;
    case 'mexFullrate'
        subbandFlag = false;
        mexFlag = true;
    case 'mexSubband'
        subbandFlag = true;
        parallelFlag = false;
        allChannelsFlag = false;
        mexFlag = true;
    case 'mexSubbandAllChannels'
        subbandFlag = true;
        parallelFlag = false;
        allChannelsFlag = true;
        mexFlag = true;
    case 'mexSubbandParallel'
        subbandFlag = true;
        parallelFlag = true;
        allChannelsFlag = false;
        mexFlag = true;
    otherwise
        error('Unknown algorithm specified');
end

%% Check combination of input parameters
switch task
    case 'speedTest'
        if ~mexFlag
            error('The speed tests expect to be run using the compiled mex files.');
        end
    case 'plotMisalignment'
        if mexFlag
            error('The misalignment analysis is only available using the MATLAB implementation.');
        end
    case 'plotDirectivity'
        if strcmp(evaluationEnvironment,'room')
            error('No directivity plot option is available for the room measurements.');
        end
    case 'plotTimeContrast'
    case 'plotPressureSpectra'

    case 'calculatePeaq'
    otherwise
        error('Selected task is not an option')
end

%% Preproccess input signal(s)
if musicFlag
    if ~exist('MiscHelperFunctions\music.wav','file')
        error('No music signal supplied. Please place a music signal in MiscHelperFunctions/music.wav')
    end
    [inSigA,fs] = audioread('MiscHelperFunctions/music.wav');
    inSigA = inSigA(31*fs:end);
    sigLen = 30*fs;
    inSigA = inSigA(1:sigLen);
    inSigB = inSigA*0;  
else
    fs = 48e3;
    sigLen = 30*fs;
    rng(1000)
    inSigA = randn(sigLen + fs,1)/10;
    inSigB = zeros(sigLen + fs,1);
    sigLen = 30*fs;
end

% Create x-over network
lowFreq{1} = 100;
lowFreq{2} = 1500;
lowFreq{3} = 4000;
highFreq{1} = 1500;
highFreq{2} = 4000;
highFreq{3} = 10000;
load('MiscHelperFunctions\crossoverNetwork.mat','sosx');
if ~strcmp(evaluationEnvironment, 'simulation')
    load('MiscHelperFunctions\parametricEqualization.mat','Bw','Aw','Bm1','Am1','Bm2','Am2','Bt','At');
    for i = 1:numberOfCrossoverBands
        inputSignalA{i} = sosfilt(sosx{i+1},inSigA);
        inputSignalB{i} = sosfilt(sosx{i+1},inSigB);
        if i == 1
            wGain = 10^(-1/20);
%           Apply biquad with center frequency 1500 Hz
%           Gain: 4 dB
%           Q: 1
%           Type: High shelving filter
            inputSignalA{i} = filter(Bw,Aw,inputSignalA{i}) * wGain;
            inputSignalB{i} = filter(Bw,Aw,inputSignalB{i}) * wGain;
        elseif i == 2
%           Apply biquad with center frequency 2700 Hz
%           Gain: -5 dB
%           Quality factor: 2
%           Type: Peak/Dip
            inputSignalA{i} = filter(Bm1,Am1,inputSignalA{i});
            inputSignalB{i} = filter(Bm1,Am1,inputSignalB{i});
%           Apply biquad with center frequency 2200 Hz
%           Gain: 4 dB
%           Quality factor: 2
%           Type: Peak/Dip
            inputSignalA{i} = filter(Bm2,Am2,inputSignalA{i});
            inputSignalB{i} = filter(Bm2,Am2,inputSignalB{i});
        else
            tGain = 10^(6/20);
%           Apply biquad with center frequency 2500 Hz
%           Gain: -6 dB
%           Quality factor: 1/sqrt(2)
%           Type: Low shelving filter
            inputSignalA{i} = filter(Bt,At,inputSignalA{i}) * tGain;
            inputSignalB{i} = filter(Bt,At,inputSignalB{i}) * tGain;
        end
    end
else
    for i = 1:numberOfCrossoverBands
        inputSignalA{i} = sosfilt(sosx{i+1},inSigA);
        inputSignalB{i} = sosfilt(sosx{i+1},inSigB);
    end
end

%% Adjust which subband channels should be active
if subbandFlag
    subbandCenterFrequencies = zeros(numberOfChannels,1);
    for i = 0:numberOfChannels-1
        subbandCenterFrequencies(i+1) = fs/numberOfChannels*(i+0.5);
    end
    subbandLowerCutoffFrequencies = subbandCenterFrequencies-0.5*fs/numberOfChannels;
    subbandLowerCutoffFrequencies(subbandLowerCutoffFrequencies<0) = 0;
    subbandHigherCutoffFrequencies = subbandCenterFrequencies+0.5*fs/numberOfChannels;
    
    oct = 1/3;
    for band = 1:numberOfCrossoverBands
        if allChannelsFlag
            gainFactors{band} = ones(numberOfChannels/2,1);
        else
            gainFactors{band} = zeros(numberOfChannels/2,1);
            for channel = 1:numberOfChannels/2
                if (subbandHigherCutoffFrequencies(channel) > lowFreq{band}/2^oct) && (subbandLowerCutoffFrequencies(channel) < highFreq{band}*2^oct)
                    gainFactors{band}(channel) = 1;
                end
            end
        end
    end
end
%% Select microphones
zoneAIdx = 15:17;
zoneBIdx = 21:23;

load(['Decompose rirs into subbands\DecomposedRirs\' setupEnvironment '_K=' int2str(numberOfChannels) '_D=' int2str(decimationFactor) '_Lp=' int2str(prototypeFilterLength)]);

%% Generate the subband signals
for bandIdx = 1:numberOfCrossoverBands
    sourceReferenceIdx = round(size(decomposedRirs{bandIdx},3)/2);
    numberOfLoudspeakers = size(decomposedRirs{bandIdx},3);
    firTaps = 300;
    modellingDelay = 150;
    
    if subbandFlag % DO SUBBAND FILTERING %%
        % Prepare inputs for adaptive processing
        firTaps = ceil(firTaps/decimationFactor);
        modellingDelay = ceil(modellingDelay/decimationFactor);

        disp(['Determining subband signals for step-size computation. Loudspeaker type idx = ' int2str(bandIdx) '/3.'])
        if mexFlag
            varStr1 = ['testAnalysisFilterbankSpeed_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            paramStr = ['(prototypeFilter, inputSignalA{' int2str(bandIdx) '});'];
            analyzedSignals = eval([varStr1 paramStr]);        
        else
            tmpAnalysisFilterbank = analysisFilterbankFast(prototypeFilter, numberOfChannels, decimationFactor);
            numberOfBlocks = floor(sigLen/decimationFactor);
            analyzedSignals = zeros(numberOfChannels/2, numberOfBlocks);
            for i = 1:numberOfBlocks
                idx = (i-1)*decimationFactor + (1:decimationFactor);
                tmpSamples = tmpAnalysisFilterbank.processInputBuffer(inputSignalA{bandIdx}(idx));
                analyzedSignals(:,i) = tmpSamples;
            end
        end
        subbandSignalLength = size(analyzedSignals,2);
        disp('Subband signals determined.')

        firstTimeThrough = true;
        for cIdx = 1:numberOfChannels/2
            zerosSignal = zeros(size(analyzedSignals,2),1);
            if(gainFactors{bandIdx}(cIdx))
                disp(['Determining max step-size for subband idx = ' int2str(cIdx) '/' num2str(numberOfChannels/2) ' | loudspeaker type idx = ' int2str(bandIdx) '/3'])
                [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon,  wiener1, wiener2] = prepareInputsForLeakyNlms(analyzedSignals(cIdx,:).', zerosSignal, decomposedRirs{bandIdx}(:,:,:,cIdx), firTaps, sourceReferenceIdx, zoneAIdx, zoneBIdx, modellingDelay);
                disp(['Max step-size determined.'])
                if firstTimeThrough
                    firstTimeThrough = false;
                    % Initialize structures for calculating the subband
                    % adaptive filtering
                    brightRirArray = zeros(size(brightRir,1),size(brightRir,2), size(brightRir,3),numberOfChannels/2);
                    darkRirArray = brightRirArray;
                    targetFilterSpectraArray = zeros(size(brightRir,1),size(brightRir,2), numberOfChannels/2);
                    stepRangeArray = zeros(2,numberOfChannels/2);
                    regParameterArray = zeros(numberOfChannels/2,1);
                    signalEpsilonArray = zeros(numberOfChannels/2,1);
                    % Initialize structures for evaluating the tracking
                    % performance
                    if ~mexFlag
                        hopSize = fftSize/2;
                        numberOfBlocks = floor(size(analyzedSignals,2)/hopSize);
                        numberOfSubbandSamples = numberOfBlocks*hopSize;
                        misalignment{bandIdx} = zeros(numberOfBlocks, numberOfChannels/2);
                        loudspeakerSubbandSignal = zeros(numberOfChannels/2, subbandSignalLength, numberOfLoudspeakers);
                        wienerList1 = zeros(fftSize, numberOfLoudspeakers, numberOfChannels/2,1);
                        wienerList2 = zeros(fftSize, numberOfLoudspeakers, numberOfChannels/2,1);
                    end
                end
                brightRirArray(:,:,:,cIdx) = brightRir;
                darkRirArray(:,:,:,cIdx) = darkRir;
                targetFilterSpectraArray(:,:,cIdx) = targetFilterSpectra;
                stepRangeArray(:,cIdx) = stepRange;
                regParameterArray(cIdx) = regParameter;
                signalEpsilonArray(cIdx) = signalEpsilon;
                if ~mexFlag
                    wienerList1(:,:, cIdx) = wiener1;
                    wienerList2(:,:, cIdx) = wiener2;
                end
            end
        end
        
        if mexFlag % DO SUBBAND MEX PROCESSING %%
            if parallelFlag
                varStr1 = ['testSubbandAdaptiveFilteringMicsParallel_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            else
                varStr1 = ['testSubbandAdaptiveFilteringMultipleMics_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            end
            paramStr = ['(prototypeFilter, inputSignalA{' int2str(bandIdx) '}, brightRirArray, darkRirArray, targetFilterSpectraArray, signalEpsilonArray, regParameterArray, stepRangeArray, logical(gainFactors{' int2str(bandIdx) '}));'];
        
            switch bandIdx
                case 1
                    varStr2 = ['_src9'];
                case 2
                    varStr2 = '_src11';
                case 3
                    varStr2 = '_src7';
            end
            for repeatIdx = 1:numberOfRepeats
                disp(['Performing adaptive processing for loudspeaker type idx = ' int2str(bandIdx) '/3.'])
                [fullMexLoudspeakerSignal, durationInNanoseconds] = eval([varStr1 varStr2 paramStr]);
                disp(['Mex computations took ' num2str(durationInNanoseconds*1e-9,'%.2f') ' s'])
                duration{bandIdx}(repeatIdx) = durationInNanoseconds*1e-9;
            end
            tmpLoudspeakerSignals{bandIdx} = fullMexLoudspeakerSignal;
        else %% DO SUBBAND MATLAB PROCESSING %%
            % Adaptive processing
            for cIdx = 1:numberOfChannels/2
                tmpWiener = wienerList1(:,:, cIdx);
                tmpAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, firTaps, brightRirArray(:,:,:,cIdx), darkRirArray(:,:,:,cIdx), targetFilterSpectraArray(:,:,cIdx), stepRangeArray(:,cIdx), maxIterations, regParameterArray(cIdx), signalEpsilonArray(cIdx));
                analysisSignal = analyzedSignals(cIdx,:).';
                if gainFactors{bandIdx}(cIdx) 
                    disp(['Performing adaptive processing for subband idx = ' int2str(cIdx) '/' int2str(numberOfChannels/2) ' | loudspeaker type idx = ' int2str(bandIdx) '/3.'])
                    tic
                    for i = 1:numberOfBlocks
                        if i == floor(numberOfBlocks/2)
                            if switchZonesFlag
                                tmpAdaptiveBeamformer.switchZones();
                                tmpWiener = wienerList2(:,:, cIdx);
                            end
                        end
                        idx = (i-1)*hopSize + (1:hopSize);
                        tmpSamples = tmpAdaptiveBeamformer.processInputBuffer(analysisSignal(idx));
                        loudspeakerSubbandSignal(cIdx,idx,:) = tmpSamples;
                        misalignment{bandIdx}(i, cIdx) = norm(tmpWiener-tmpAdaptiveBeamformer.getFilterSpectra(),'fro')/norm(tmpWiener,'fro');
                    end
                    disp(['Matlab OO computations took ' num2str(toc,'%.2f') ' s'])
                end
                if cIdx == numberOfChannels/2
                    numberOfBlocks = size(loudspeakerSubbandSignal,2);
                    ooLoudspeakerSignal = zeros(numberOfBlocks*decimationFactor,numberOfLoudspeakers);
                    for lIdx = 1:numberOfLoudspeakers
                        ooLoudspeakerSignal(:,lIdx) = synthesizeSignal(loudspeakerSubbandSignal(:,:,lIdx), prototypeFilter, numberOfChannels, decimationFactor);
                    end
                    tmpLoudspeakerSignals{bandIdx} = ooLoudspeakerSignal;
                end
            end
        end

    else %% DO FULLRATE PROCESSING %%
        sourceReferenceIdx = round(size(decomposedRirs{bandIdx},3)/2);
        numberOfLoudspeakers = size(decomposedRirs{bandIdx},3);
        firTaps = 300;
        modellingDelay = 150;
        
        % Prepare inputs for the adaptive processing
        zerosSignal = zeros(length(inputSignalA{bandIdx}),1);
        disp(['Determining max step-size for loudspeaker type idx = ' int2str(bandIdx) '/3'])
        [fftSize, brightRir, darkRir, targetFilterSpectra, stepRange, regParameter, signalEpsilon,  wiener1, wiener2] = prepareInputsForLeakyNlms(inputSignalA{bandIdx}, zerosSignal, IRs{bandIdx}(:,:,:), firTaps, sourceReferenceIdx, zoneAIdx, zoneBIdx, modellingDelay);
        disp(['Max step-size determined.'])
        hopSize = fftSize/2;
        numberOfBlocks = floor(sigLen/hopSize);

        if mexFlag %% DO FULLRATE MEX PROCESSING %%
            varStr1 = ['testAdaptiveFilteringMultiMicsSpeedReal'];
            paramStr = ['(inputSignalA{' int2str(bandIdx) '}, brightRir, darkRir, targetFilterSpectra, signalEpsilon, regParameter, stepRange);'];
           
            switch bandIdx
                case 1
                    varStr2 = ['_src9'];
                case 2
                    varStr2 = '_src11';
                case 3
                    varStr2 = '_src7';
            end
            for repeatIdx = 1:numberOfRepeats
                disp(['Performing adaptive processing for loudspeaker type idx = ' int2str(bandIdx) '/3.'])
                [fullMexLoudspeakerSignal, durationInNanoseconds] = eval([varStr1 varStr2 paramStr]);
                disp(['Mex computations took ' num2str(durationInNanoseconds*1e-9,'%.2f') ' s'])
                duration{bandIdx}(repeatIdx) = durationInNanoseconds*1e-9;
                tmpLoudspeakerSignals{bandIdx} = fullMexLoudspeakerSignal;
            end
        else %% DO FULLRATE MATLAB PROCESSING %%
            tmpAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, firTaps, brightRir, darkRir, targetFilterSpectra, stepRange, maxIterations, regParameter, signalEpsilon);
            loudspeakerSignal = zeros(numberOfBlocks*hopSize, numberOfLoudspeakers);
            misalignment{bandIdx} = zeros(numberOfBlocks,1);
            tmpWiener = wiener1;
            analysisSignal = inputSignalA{bandIdx};
            disp(['Performing adaptive processing for loudspeaker type idx = ' int2str(bandIdx) '/3.'])
            tic
            for i = 1:numberOfBlocks
                if i == floor(numberOfBlocks/2)
                    if switchZonesFlag
                        tmpAdaptiveBeamformer.switchZones();
                        tmpWiener = wiener2;
                    end
                end
                idx = (i-1)*hopSize + (1:hopSize);
                tmpSamples = tmpAdaptiveBeamformer.processInputBuffer(analysisSignal(idx));
                loudspeakerSignal(idx,:) = tmpSamples;
                misalignment{bandIdx}(i) = norm(tmpWiener-tmpAdaptiveBeamformer.getFilterSpectra(),'fro')/norm(tmpWiener,'fro');
            end
            disp(['Matlab OO computations took ' num2str(toc,'%.2f') ' s'])
            tmpLoudspeakerSignals{bandIdx} = loudspeakerSignal;
        end
    end
end

%% Post-process results
switch task
    case 'speedTest'
        saveFolderName = 'AdaptiveFilteringTests/TimeResults/';
        if subbandFlag
            if parallelFlag
                fname = ['processingDurationParallel_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            elseif allChannelsFlag
                fname = ['processingDurationAllChannels_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            else
                fname = ['processingDuration_K' int2str(numberOfChannels) '_D' int2str(decimationFactor) '_Lp' int2str(prototypeFilterLength)];
            end
        else
            fname = ['processingDurationReal'];
        end
        save([saveFolderName fname], "duration");

    case 'plotMisalignment'
        if mexFlag
            error('The analysis for the misalignment is only available for the matlab implementation.');
        end
        timeIdx = ((1:length(misalignment{1}))-1)*hopSize*decimationFactor/fs;
        lineType{1} = '-';
        lineType{2} = '--';
        lineType{3} = ':';
        figure
        for bandIdx = 1:numberOfCrossoverBands
            if subbandFlag
                for cIdx = 1:numberOfChannels/2
                    if gainFactors{bandIdx}(cIdx)
                        plot(timeIdx, 20*log10(misalignment{bandIdx}(:,cIdx)), lineType{bandIdx},'color','blue');
                        hold on; grid on
                    end
                end
            else
                plot(timeIdx, 20*log10(misalignment{bandIdx}), lineType{bandIdx},'color','blue');
                hold on; grid on
            end
        end
        xlabel('Time [s]'); ylabel('NMSE [dB]')

    case 'plotDirectivity'
        switch evaluationEnvironment
            case 'simulation'
            case 'anechoic'
                load('TransferFunctions/Prototype2_Setup1_FreeFieldResponses')
                tmpIRs = resp;
                midIdx = [21:29];
                subTwIdx = [9:19];
                twIdx = [1:7];
                IRs{1} = tmpIRs(:,:,midIdx);
                IRs{2} = tmpIRs(:,:,subTwIdx);
                IRs{3} = tmpIRs(:,:,twIdx);
            otherwise
                error(['Unknow EvaluationEnvironment specified for the task ' task]);
        end
        PlotDirectivityResults(IRs, tmpLoudspeakerSignals)

    case 'plotTimeContrast'
        switch evaluationEnvironment
            case 'simulation'
            case 'anechoic'
                load('TransferFunctions/Prototype2_Setup1_FreeFieldResponses')
                tmpIRs = resp;
                midIdx = [21:29];
                subTwIdx = [9:19];
                twIdx = [1:7];
                IRs{1} = tmpIRs(:,:,midIdx);
                IRs{2} = tmpIRs(:,:,subTwIdx);
                IRs{3} = tmpIRs(:,:,twIdx);
            case 'room'
                load('TransferFunctions/roomEvaluationMeasurementsZoneA.mat');
                load('TransferFunctions/roomEvaluationMeasurementsZoneB.mat');
                tmpIRs = cat(2,EvaluationRirZoneA,EvaluationRirZoneB);
                midIdx = [21:29];
                subTwIdx = [9:19];
                twIdx = [1:7];
                IRs{1} = tmpIRs(:,:,midIdx);
                IRs{2} = tmpIRs(:,:,subTwIdx);
                IRs{3} = tmpIRs(:,:,twIdx);
                zoneAIdx = 1:12;
                zoneBIdx = 13:24;
            otherwise
                error(['Unknow EvaluationEnvironment specified for the task ' task]);
        end
        PlotTimeContrast(IRs, tmpLoudspeakerSignals, zoneAIdx, zoneBIdx);

    case 'plotPressureSpectra'
        switch evaluationEnvironment
            case 'simulation'
            case 'anechoic'
                load('TransferFunctions/Prototype2_Setup1_FreeFieldResponses')
                tmpIRs = resp;
                midIdx = [21:29];
                subTwIdx = [9:19];
                twIdx = [1:7];
                IRs{1} = tmpIRs(:,:,midIdx);
                IRs{2} = tmpIRs(:,:,subTwIdx);
                IRs{3} = tmpIRs(:,:,twIdx);
            case 'room'
                load('TransferFunctions/roomEvaluationMeasurementsZoneA.mat');
                load('TransferFunctions/roomEvaluationMeasurementsZoneB.mat');
                tmpIRs = cat(2,EvaluationRirZoneA,EvaluationRirZoneB);
                midIdx = [21:29];
                subTwIdx = [9:19];
                twIdx = [1:7];
                IRs{1} = tmpIRs(:,:,midIdx);
                IRs{2} = tmpIRs(:,:,subTwIdx);
                IRs{3} = tmpIRs(:,:,twIdx);
                zoneAIdx = 1:12;
                zoneBIdx = 13:24;
            otherwise
                error(['Unknow EvaluationEnvironment specified for the task ' task]);
        end
        PlotPressureSpectra(IRs, tmpLoudspeakerSignals, zoneAIdx, zoneBIdx);

    case 'calculatePeaq'
        ReferenceSignal = zeros(length(inputSignalA{1}),1);
        for idx = 1:length(inputSignalA)
            ReferenceSignal = ReferenceSignal + inputSignalA{idx};
        end
        Pressure = PredictPressure(tmpLoudspeakerSignals, IRs, round(mean(zoneAIdx)));        
        Latency = varargin{1};
        
        ReferenceSignal = ReferenceSignal/sqrt(mean(ReferenceSignal.^2));
        Pressure = Pressure(Latency:end-fs)/sqrt(mean(Pressure(Latency:end-fs).^2));
        OutputSignalLength = length(Pressure);
        ReferenceSignal = ReferenceSignal(1:OutputSignalLength);
        Scale = max(abs(ReferenceSignal)) * 1.1;
        ReferenceSignal = ReferenceSignal/Scale;
        Pressure = Pressure / Scale;
        
        audiowrite('ref.wav',ReferenceSignal(3*fs:end),fs);
        audiowrite('test.wav',Pressure(3*fs:end),fs);
        
        [PeaqScore, ~] = PQevalAudio_fn('ref.wav', 'test.wav');
        varargout{1} = PeaqScore;
    otherwise
        error('Specified task not recognized as a valid post-processing option')
end


end
