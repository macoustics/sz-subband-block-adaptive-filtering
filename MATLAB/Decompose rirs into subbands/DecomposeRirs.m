function [SubbandRir] = DecomposeRirs(AnalysisFilters, DecimationFactor, Rir)
% Function for estimating the subband decomposition of room impulse
% responses.
%
% =============================
% INPUTS:
% -----------------------------
% AnalysisFilters   2D array
% size = [PrototypeFilterLength, NumberOfChannels]
% DecimationFactor  int
% Rir               cell array   Room impulse responses
% size(Rir) = [NumberOfCrossoverBands,1]
% size(Rir{x}) = [RirLength, NumberOfMicrophones, NumberOfSources]
%
% =============================
% OUTPUTS:
% -----------------------------
% DecomposedRir     cell array  Room impulse responses decomposed into
%                               subband systems.
% size(DecomposedRir) = [NumberOfCrossoverBands,1]
% size(DecomposedRir{x}) = [SubsystemResponseLength, NumberOfMicrophones, 
% NumberOfSources, NumberOfChannels].

%% Analysis
PrototypeFilterLength = size(AnalysisFilters,1);
NumberOfChannels = size(AnalysisFilters,2);
NumberOfCrossoverBands = length(Rir);
for band = 1:NumberOfCrossoverBands
    ss = size(Rir{band});
    RirLength = ss(1);
    NumberOfMicrophones = ss(2);
    NumberOfSources = ss(3);
    SubsampledConvolutionLength = ceil((RirLength + PrototypeFilterLength-1)/DecimationFactor); 
    SubsampledAnalysisFilterLength = ceil(PrototypeFilterLength/DecimationFactor);
    SubsystemResponseLength = SubsampledConvolutionLength - ceil(PrototypeFilterLength/DecimationFactor) + 1;
    AnalyzedRir = zeros(SubsystemResponseLength,NumberOfMicrophones,NumberOfSources,NumberOfChannels);
    for src = 1:NumberOfSources
        for mic = 1:NumberOfMicrophones
            disp(['Band = ' int2str(band) '/' int2str(NumberOfCrossoverBands) ' | Src = ' int2str(src) '/' int2str(NumberOfSources) ' | Mic = ' int2str(mic) '/' int2str(NumberOfMicrophones)])
            tmpIR = Rir{band}(:,mic,src);
            tmpIR = tmpIR(:);
            for bank = 1:NumberOfChannels
                ak = downsample(AnalysisFilters(:,bank),DecimationFactor);
                Uk = toeplitz([ak; zeros(SubsampledConvolutionLength-SubsampledAnalysisFilterLength,1)], [ak(1), zeros(1,SubsystemResponseLength-1)]);
                TargetSubsystemResponse = downsample(conv(tmpIR,AnalysisFilters(:,bank),'full'), DecimationFactor);
                tmp = Uk'*Uk;
                if cond(tmp) > 1e5
                    AnalyzedRir(:,mic,src,bank) = (tmp + 1e-5*eye(length(tmp))*norm(tmp))\(Uk'*TargetSubsystemResponse);
                else
                    AnalyzedRir(:,mic,src,bank) = Uk\TargetSubsystemResponse;
                end
            end
        end
    end
    SubbandRir{band} = AnalyzedRir;
end

