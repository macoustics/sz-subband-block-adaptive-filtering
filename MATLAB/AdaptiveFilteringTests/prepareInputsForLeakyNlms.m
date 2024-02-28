function [fftSize, brightRir, darkRir, targetFilter, stepRange, regParameter, signalEpsilon, Wiener1, Wiener2] = prepareInputsForLeakyNlms(SignalA, SignalB, RoomImpulseResponses, FirTaps, SourceReferenceIdx, MicReferenceIdxA, MicReferenceIdxB, ModellingDelay)
% INPUTS
% --------------------------
% SignalA             vector      Audio signal played in zone A
% SignalB             vector      Audio signal played in zone B
% RoomImpulseResponsess   
%                     Array       Room Impulse responses
%                                 size = [nSamples, nMiccs, nSrcs]
% fs                  Int         Sampling frequency [Hz]
% FirTaps             Int         Length of control filters
% SourceReferenceIdx  Int         Index to reference loudspeaker loudspeaker 
% MicReferenceIdxA    Int         Index to reference microphone for zone A
% MicReferenceIdxB    Int         Index to reference microphone for zone B
% TestLength          Int         Length of input signal to be processed
% gamma               Scalar      Forgetting factor for recursive estimates
%
% OUTPUTS
% --------------------------
% brightRir   
%                     Array       Room Impulse responses for bright zone
%                                 size = [nSamples, nSrcs, nBrightMics]
% darkRir   
%                     Array       Room Impulse responses for dark zone
%                                 size = [nSamples, nSrcs, nDarkMics]
% targetFilter        Array       FIR filters describing the target
%                                 pressure in the bright zone.
%                                 size = [nSamples, nSrcs]
% stepRange           Vector      [minStepSize; maxStepSize]
% regParameter        Scalar      Regularization parameter calculated from
%                                 the particular room impulse responses.
% signalEpsilon       Scalar      Regularization parameter dependent on the
%                                 particular input signal.


%% Tuning parameters

kRegularizationScaling = 1e-4;%1e-4;



%% Initialize parameters

% Filter lengths (currently we assume same lengths for all filters)
kNumberOfSources = size(RoomImpulseResponses,3); % Number of channels (loud speakers)
kNumberOfMics = length(MicReferenceIdxA); % Number of microphones in each zone
    
kRoomImpulseResponseLength = size(RoomImpulseResponses,1); % Length of RIRs for both bright and dark zone
kConvolutionLength = kRoomImpulseResponseLength + FirTaps - 1;
kNfft = 2*(kConvolutionLength) - 1; % The fft should include both control filter and acoustic impulse response
kNfft = 2^nextpow2(kNfft);
fftSize = kNfft;
kHopSize = kNfft/2; % The output only takes the control filter into account
disp(['FFT block size = ' int2str(kNfft)]);

% Introduce a parallel filter structure for determining a target sound
% field in the bright zone, which we can use as a least-squares term in the
% cost function.
mTargetFilter = zeros(FirTaps,kNumberOfSources); % Equalization filters

% Make the equalization filters into modelling delays, which doesn't
% change with time. The modelling delay allows the filters to have
% pre-ringing to align the phase of each component.
mTargetFilter(ModellingDelay, SourceReferenceIdx) = 1;
targetFilter = fft(mTargetFilter, fftSize, 1);



% Matrices containing RIRs for the bright and dark zones.
% Row k, is the filter for the k-th time slot.
% if lowIRs
% mRoomImpulseResponses =  permute(RoomImpulseResponses, [2,1,3]);
mRoomImpulseResponses1 = RoomImpulseResponses(1:kRoomImpulseResponseLength,MicReferenceIdxA,:);
mRoomImpulseResponses2 = RoomImpulseResponses(1:kRoomImpulseResponseLength,MicReferenceIdxB,:);
% Permute from [Samples, Mics, Srcs] to [Samples, Srcs, Mics]
mRoomImpulseResponses1 = permute(mRoomImpulseResponses1, [1,3,2]);
mRoomImpulseResponses2 = permute(mRoomImpulseResponses2, [1,3,2]);

% FFT to required blockSize
mRoomTransferFunctions1 = fft(mRoomImpulseResponses1,kNfft,1);
mRoomTransferFunctions2 = fft(mRoomImpulseResponses2,kNfft,1);

mTest = zeros(kNfft,1);
for idx = 1:kNfft
    if kNumberOfMics > 1
        mTest(idx) = norm(cat(1,squeeze(mRoomTransferFunctions1(idx,:,:)).',squeeze(mRoomTransferFunctions2(idx,:,:)).'));
    else
        mTest(idx) = norm(cat(1,squeeze(mRoomTransferFunctions1(idx,:,:)),squeeze(mRoomTransferFunctions2(idx,:,:))));
    end
end
% normConstant = 0;
% for fIdx = 1:kNfft
%     HH = zeros(kNumberOfSources);
%     for m=1:kNumberOfMics
%         tf = mRoomTransferFunctions1(fIdx,:,m);
%         tf = tf(:);
%         HH = HH + conj(tf)*(tf.');
%         tf = mRoomTransferFunctions2(fIdx,:,m);
%         tf = tf(:);
%         HH = HH + conj(tf)*(tf.');
%     end
%     normConstantTmp = norm(HH);
%     if normConstantTmp > normConstant
%         normConstant = normConstantTmp;
%     end
% end
% clear HH tf

disp(['Before normalization: normRTF1 = ' num2str(norm(mRoomTransferFunctions1(:)), '%.2e') ' normRTF2 = ' num2str(norm(mRoomTransferFunctions2(:)), '%.2e')])
mRoomTransferFunctions1 = mRoomTransferFunctions1/max(mTest);
mRoomTransferFunctions2 = mRoomTransferFunctions2/max(mTest);
disp(['After normalization: normRTF1 = ' num2str(norm(mRoomTransferFunctions1(:)), '%.2e') ' normRTF2 = ' num2str(norm(mRoomTransferFunctions2(:)), '%.2e')])
brightRir = mRoomTransferFunctions1;
darkRir = mRoomTransferFunctions2;


%% Pre-calculate matrices
kRegularizationConstant = 0;
for fIdx = 1:kNfft
    HH = zeros(kNumberOfSources);
    for m=1:kNumberOfMics
        tf = mRoomTransferFunctions1(fIdx,:,m);
        tf = tf(:);
        HH = HH + conj(tf)*(tf.');
        tf = mRoomTransferFunctions2(fIdx,:,m);
        tf = tf(:);
        HH = HH + conj(tf)*(tf.');
    end
    RegularizationConstantTmp = kRegularizationScaling*norm(HH);
    if RegularizationConstantTmp > kRegularizationConstant
        kRegularizationConstant = RegularizationConstantTmp;
    end
end

% keyboard
kRegularizationConstant = kRegularizationConstant;
regParameter = kRegularizationConstant;


%% Determine max step length
% size(mRoomTransferFunctions) = [kNfft,NumberOfSources,NumberOfMics]
% keyboard
% Est. signal cross-correlation matrix
kSamplePosition = 0;
AvgNumber = 0;
W = dftmtx(kNfft);
G01 = 1/2*W*[zeros(kNfft/2), zeros(kNfft/2); zeros(kNfft/2), eye(kNfft/2)]*W'/kNfft;

SignalCovariance = zeros(kNfft);
while kSamplePosition + kNfft <= length(SignalA)
    % Get current signal vector
    mInputSignal1 = SignalA(kSamplePosition + (1:kNfft));
    mInputSpectrum1 = fft(mInputSignal1,kNfft);
    SignalCovariance = SignalCovariance + conj(mInputSpectrum1)*mInputSpectrum1.';
    kSamplePosition = kSamplePosition + kHopSize;
    AvgNumber = AvgNumber + 1;
end
SignalCovariance = SignalCovariance/AvgNumber;
signalEpsilon = 1e-2*max(max(abs(diag(SignalCovariance))));
invLambda = diag(1./(diag(SignalCovariance) + signalEpsilon));
RfG = (SignalCovariance .* G01) * invLambda;
R = zeros(2*kNumberOfMics*kNfft,kNumberOfSources*kNfft);
for sIdx = 1:kNumberOfSources
    for mIdx = 1:kNumberOfMics
        idx = (sIdx-1)*kNfft + (1:kNfft);
        R((mIdx-1)*kNfft + (1:kNfft),idx) = diag(mRoomTransferFunctions1(:,sIdx,mIdx));
        R(kNumberOfMics*kNfft+(mIdx-1)*kNfft+(1:kNfft),idx) = diag(mRoomTransferFunctions2(:,sIdx,mIdx));
    end
end

I = sparse(eye(kNumberOfMics*2));
RfMult = kron(I,RfG);
B = R'*RfMult*R;

% keyboard
A = B + kRegularizationConstant*eye(length(B));
% keyboard
% condA = cond(A);
Z10 = [eye(FirTaps); zeros(kNfft-FirTaps, FirTaps)];
WZ10 = W*Z10;
Rc = kron(eye(kNumberOfSources),WZ10')*A'*kron(eye(kNumberOfSources),WZ10)/kNfft;
tmp3 = eigs(Rc,1);

kMaxStepLength = 2/abs(tmp3)/2;
kMinStepLength = kMaxStepLength/1e3;

stepRange = [kMinStepLength, kMaxStepLength];

%% Determine Wiener solution
G10 = W * [eye(FirTaps); zeros(kNfft-FirTaps,FirTaps)];
G10 = kron(eye(kNumberOfSources),G10);

Hessian = zeros(kNfft*kNumberOfSources);
RHS1 = zeros(kNfft*kNumberOfSources,1);
RHS2 = zeros(kNfft*kNumberOfSources,1);
mTargetFilterSpectrum = targetFilter(:,SourceReferenceIdx);
for mIdx = 1:kNumberOfMics
     rTmp = zeros(kNfft,kNfft*kNumberOfSources);
     for sIdx = 1:kNumberOfSources
         idx = (sIdx-1)*kNfft + (1:kNfft);
         rTmp(:,idx) = diag(mRoomTransferFunctions1(:,sIdx,mIdx));
     end
     Hessian = Hessian + rTmp'*RfG*rTmp;
     RHS1 = RHS1 + rTmp'*RfG*(mRoomTransferFunctions1(:,SourceReferenceIdx,mIdx).*mTargetFilterSpectrum);
     rTmp = zeros(kNfft,kNfft*kNumberOfSources);
     for sIdx = 1:kNumberOfSources
         idx = (sIdx-1)*kNfft + (1:kNfft);
         rTmp(:,idx) = diag(mRoomTransferFunctions2(:,sIdx,mIdx));
     end
     Hessian = Hessian + rTmp'*RfG*rTmp;
     RHS2 = RHS2 + rTmp'*RfG*(mRoomTransferFunctions2(:,SourceReferenceIdx,mIdx).*mTargetFilterSpectrum);
end
% Hessian = Hessian + 1/2*kRegularizationConstant*kron(eye(kNumberOfSources), diag(diag(SignalCovariance)));
% Wwiener = Hessian\RHS;
% Wwiener = TruncateFilters(Wwiener, kNumberOfSources, FirTaps);

% Hessian = Hessian + 1/2*kRegularizationConstant*kron(eye(kNumberOfSources), diag(diag(SignalCovariance)));
Hessian = Hessian + kRegularizationConstant*eye(size(Hessian))*1;
Hessian = G10'*Hessian*G10;
RHS1 = G10'*RHS1;
Wiener1 = Hessian\RHS1;
Wiener1 = reshape(Wiener1,FirTaps,kNumberOfSources);
Wiener1 = fft(Wiener1, kNfft, 1);
% keyboard

RHS2 = G10'*RHS2;
Wiener2 = Hessian\RHS2;
Wiener2 = reshape(Wiener2,FirTaps,kNumberOfSources);
Wiener2 = fft(Wiener2, kNfft, 1);

disp(['kMaxStepLength = ' num2str(kMaxStepLength,'%.2e') ' RegConst = ' num2str(kRegularizationConstant,'%.2e') ' signalEpsilon = ' num2str(signalEpsilon,'%.2e')])








