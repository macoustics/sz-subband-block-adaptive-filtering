classdef leakyNlmsAdaptiveBeamformer < handle
    properties
        %% General settings
        m_numberOfBrightMicrophones = 0;
        m_numberOfDarkMicrophones = 0;
        m_numberOfLoudspeakers = 0;
        m_filterLength = 0;
        m_fftSize = 0;
        m_brightRir = [];
        m_targetFilterSpectrum = [];
        m_darkRir = [];
        m_filterSpectrum = [];
        m_powerSpectrum = [];
        m_inputSignal = [];
        m_inputSpectrum = [];
        m_regularizationFactor = 0;
        m_signalEpsilon = 0;

        %% Line search settings
        m_gammaSpectrum = 1-1e-2;
        m_maxIterations = 0;
        m_stepSize = 0;
        m_minStepSize = 0;
        m_maxStepSize = 0;
        m_armijoConstant = 1e-4;
        m_stepScale = 0.1;
    end
    methods 
        function leakyNlmsAdaptiveBeamformer = leakyNlmsAdaptiveBeamformer(fftSize, filterLength, brightRir, darkRir, targetBrightFilterSpectra, stepRange, maxIterations, regParameter, signalEpsilon)
            %% Initialize the beamformer
            leakyNlmsAdaptiveBeamformer.m_fftSize = fftSize;
            leakyNlmsAdaptiveBeamformer.m_filterLength = filterLength;

            leakyNlmsAdaptiveBeamformer.m_brightRir = brightRir;
            leakyNlmsAdaptiveBeamformer.m_darkRir = darkRir;
            leakyNlmsAdaptiveBeamformer.m_targetFilterSpectrum = targetBrightFilterSpectra;
            leakyNlmsAdaptiveBeamformer.m_numberOfBrightMicrophones = size(brightRir,3);
            leakyNlmsAdaptiveBeamformer.m_numberOfDarkMicrophones = size(darkRir,3);
            leakyNlmsAdaptiveBeamformer.m_numberOfLoudspeakers = size(brightRir,2);
            leakyNlmsAdaptiveBeamformer.m_inputSignal = zeros(fftSize, 1);
            leakyNlmsAdaptiveBeamformer.m_powerSpectrum = zeros(fftSize, 1);
            leakyNlmsAdaptiveBeamformer.m_filterSpectrum = zeros(fftSize, leakyNlmsAdaptiveBeamformer.m_numberOfLoudspeakers);

            leakyNlmsAdaptiveBeamformer.m_minStepSize = stepRange(1);
            leakyNlmsAdaptiveBeamformer.m_maxStepSize = stepRange(2);
            leakyNlmsAdaptiveBeamformer.m_stepSize = leakyNlmsAdaptiveBeamformer.m_maxStepSize/100;
            leakyNlmsAdaptiveBeamformer.m_maxIterations = maxIterations;
            leakyNlmsAdaptiveBeamformer.m_regularizationFactor = regParameter;
            leakyNlmsAdaptiveBeamformer.m_signalEpsilon = signalEpsilon;
        end
        function [filterSpectra] = getFilterSpectra(obj)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            arguments (Output)
                filterSpectra (:,:) double
            end
            filterSpectra = obj.m_filterSpectrum;
        end

        function [obj] = switchZones(obj)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            arguments (Output)
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            tmp = obj.m_brightRir;
            obj.m_brightRir = obj.m_darkRir;
            obj.m_darkRir = tmp;
        end

        function [obj] = updateInputSpectrum(obj, inputBuffer)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                inputBuffer (:,1) double
            end
            arguments (Output)
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            % 50 percent overlap for the overlap save
            obj.m_inputSignal = [obj.m_inputSignal(obj.m_fftSize/2+1:end); inputBuffer];
            obj.m_inputSpectrum = fft(obj.m_inputSignal);
            obj.m_powerSpectrum = obj.m_gammaSpectrum * obj.m_powerSpectrum + (1-obj.m_gammaSpectrum)*abs(obj.m_inputSpectrum).^2;
        end

        function [outputBuffer, obj] = updateAndApplyFilters(obj, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                gradient (:,:) double
            end
            arguments (Output)
                outputBuffer (:,:) double
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            outputBuffer = zeros(obj.m_fftSize/2, obj.m_numberOfLoudspeakers);
            tmp = zeros(obj.m_fftSize,1);
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                % Update the filters via gradient descent
                obj.m_filterSpectrum(:,sIdx+1) = obj.m_filterSpectrum(:,sIdx+1) - obj.m_stepSize*gradient(:,sIdx+1);
                % Apply the filters
                tmp = ifft(obj.m_inputSpectrum.*obj.m_filterSpectrum(:,sIdx+1));
                outputBuffer(:,sIdx+1) = tmp(obj.m_fftSize/2+1:end);
            end
            
        end

        function [outputBuffer, obj] = processInputBuffer(obj, inputBuffer)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                inputBuffer (:,1) double
            end
            arguments (Output)
                outputBuffer (:,:) double
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            obj.updateInputSpectrum(inputBuffer);
            gradient = obj.calculateGradient();
            obj.calculateStepSize(gradient);
            outputBuffer = obj.updateAndApplyFilters(gradient);
        end

        %% Calculate error
        function [errorSpectrum] = calculateBrightError(obj, rir, filterSpectrum)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                rir (:,:) double
                filterSpectrum (:,:) double
            end
            arguments (Output)
                errorSpectrum (:,1) double
            end
            errorSpectrum = zeros(obj.m_fftSize,1);
            % Calc error with circular aliasing
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                errorSpectrum = errorSpectrum + obj.m_inputSpectrum.*rir(:,sIdx+1).*(obj.m_targetFilterSpectrum(:,sIdx+1) - filterSpectrum(:,sIdx+1));
            end
            % Enforce linear convolution
            errorSequence = ifft(errorSpectrum);
            errorSequence((0:obj.m_fftSize/2-1) + 1) = 0;
            errorSpectrum = 1/2*fft(errorSequence);
        end

        function [gradient] = calculateBrightErrorUpdateGradient(obj, rir, filterSpectrum, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                rir (:,:,:) double
                filterSpectrum (:,:) double
                gradient (:,:) double
            end
            arguments (Output)
                gradient (:,:) double
            end
            for mIdx = 0:obj.m_numberOfBrightMicrophones-1
                errorSpectrum = obj.calculateBrightError(rir(:,:,mIdx+1), filterSpectrum);
                % Apply error to forward model
                for sIdx = 0:obj.m_numberOfLoudspeakers-1
                    gradient(:,sIdx+1) = gradient(:,sIdx+1) - conj(obj.m_inputSpectrum.*rir(:,sIdx+1,mIdx+1)).*errorSpectrum;
                end
            end
        end

        function [errorSpectrum] = calculateDarkError(obj, rir, filterSpectrum)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                rir (:,:) double
                filterSpectrum (:,:) double
            end
            arguments (Output)
                errorSpectrum (:,1) double
            end
            errorSpectrum = zeros(obj.m_fftSize,1);
            % Calc error with circular aliasing
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                errorSpectrum = errorSpectrum + obj.m_inputSpectrum.*(- rir(:,sIdx+1).*filterSpectrum(:,sIdx+1));
            end
            % Enforce linear convolution
            errorSequence = ifft(errorSpectrum);
            errorSequence((0:obj.m_fftSize/2-1) + 1) = 0;
            errorSpectrum = 1/2*fft(errorSequence);
        end

        function [gradient] = calculateDarkErrorUpdateGradient(obj, rir, filterSpectrum, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                rir (:,:,:) double
                filterSpectrum (:,:) double
                gradient (:,:) double
            end
            arguments (Output)
                gradient (:,:) double
            end
            for mIdx = 0:obj.m_numberOfDarkMicrophones-1
                errorSpectrum = obj.calculateDarkError(rir(:,:,mIdx+1), filterSpectrum);
                % Apply error to forward model
                for sIdx = 0:obj.m_numberOfLoudspeakers-1
                    gradient(:,sIdx+1) = gradient(:,sIdx+1) - conj(obj.m_inputSpectrum.*rir(:,sIdx+1,mIdx+1)).*errorSpectrum;
                end
            end
        end

        function [gradient] = normalizeAndAddRegularization(obj, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                gradient (:,:) double
            end
            arguments (Output)
                gradient (:,:) double
            end
            % Normalize to inputPowerSpectrum
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                gradient(:,sIdx+1) = gradient(:,sIdx+1) .* 1./(obj.m_powerSpectrum + obj.m_signalEpsilon);
            end
            % Apply the "leaky" regularization factor
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                gradient(:,sIdx+1) = gradient(:,sIdx+1) + obj.m_regularizationFactor*obj.m_filterSpectrum(:,sIdx+1);
            end
        end

        function [gradient] = forceLinearConvolution(obj, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                gradient (:,:) double
            end
            arguments (Output)
                gradient (:,:) double
            end
            % Enforce linear convolution
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                tmpGradient = ifft(gradient(:,sIdx+1));
                for idx = obj.m_filterLength : obj.m_fftSize-1
                    tmpGradient(idx+1) = 0;
                end
                gradient(:,sIdx+1) = obj.m_fftSize/obj.m_filterLength*fft(tmpGradient);
            end
        end

        %% Calculate gradient
        function [gradient] = calculateGradient(obj)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
            end
            arguments (Output)
                gradient (:,:) double
            end
            gradient = zeros(obj.m_fftSize,obj.m_numberOfLoudspeakers);
            % Gradient due to the bright microphones
            [gradient] = obj.calculateBrightErrorUpdateGradient(obj.m_brightRir, obj.m_filterSpectrum, gradient);
            % Gradient due to the dark microphones 
            [gradient] = obj.calculateDarkErrorUpdateGradient(obj.m_darkRir, obj.m_filterSpectrum, gradient);
            gradient = obj.normalizeAndAddRegularization(gradient);
            gradient = obj.forceLinearConvolution(gradient);
        end

        %% Determine the step size
        function [cost] = evalCost(obj, filterSpectrum)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                filterSpectrum (:,:) double
            end
            arguments (Output)
                cost (1,1) double
            end
            cost = 0;
            % Bright zone error
            for mIdx = 0:obj.m_numberOfBrightMicrophones-1
                errorSpectrum = obj.calculateBrightError(obj.m_brightRir(:,:,mIdx+1), filterSpectrum);
                cost = cost + sum(abs(errorSpectrum).^2);
            end
            % Dark zone error
            for mIdx = 0:obj.m_numberOfDarkMicrophones-1
                errorSpectrum = obj.calculateDarkError(obj.m_darkRir(:,:,mIdx+1), filterSpectrum);
                cost = cost + sum(abs(errorSpectrum).^2);
            end
            % Filter norm (leaky term)
%             for sIdx = 0:obj.m_numberOfLoudspeakers-1
%                 tmp = ifft(obj.m_inputSpectrum .* filterSpectrum(:,sIdx+1));
%                 tmp(1:obj.m_fftSize/2) = 0;
%                 cost = cost + obj.m_regularizationFactor*obj.m_fftSize^2/4 * sum(abs(tmp).^2);
%             end
            tmp = 0;
            for sIdx = 0:obj.m_numberOfLoudspeakers-1
                tmp = tmp + 1/2*sum(abs(obj.m_filterSpectrum(:,sIdx+1)).^2);
            end
            cost = cost + obj.m_regularizationFactor*tmp;
        end

        function [stepSize] = decreasingLineSearch(obj, newCost, prevCost, tmp, stepSize, stepDirection)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                newCost (1,1) double
                prevCost (1,1) double
                tmp (1,1) double
                stepSize (1,1) double
                stepDirection (:,:) double
            end
            arguments (Output)
                stepSize (1,1) double
            end
            lineSearchItr = 0;
            while newCost > prevCost + obj.m_armijoConstant*stepSize*tmp && lineSearchItr < obj.m_maxIterations
                lineSearchItr = lineSearchItr + 1;
                stepSize = obj.m_stepScale * stepSize;
                newCost = obj.evalCost(obj.m_filterSpectrum+stepSize*stepDirection);
            end
        end

        function [stepSize] = increasingLineSearch(obj, newCost, prevCost, tmp, stepSize, stepDirection)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                newCost (1,1) double
                prevCost (1,1) double
                tmp (1,1) double
                stepSize (1,1) double
                stepDirection (:,:) double
            end
            arguments (Output)
                stepSize (1,1) double
            end
            lineSearchItr = 0;
            while newCost <= prevCost + obj.m_armijoConstant*stepSize*tmp && lineSearchItr < obj.m_maxIterations
                lineSearchItr = lineSearchItr + 1;
                stepSize = 1/obj.m_stepScale * stepSize;
                newCost = obj.evalCost(obj.m_filterSpectrum+stepSize*stepDirection);
            end
            % Choose the last stepsize satisfying the Armijo condition
            if lineSearchItr > 0
                stepSize = obj.m_stepScale*stepSize;
            end
        end

        function [stepSize] = calculateStepSize(obj, gradient)
            arguments (Input)
                obj (1,1) leakyNlmsAdaptiveBeamformer
                gradient (:,:) double
            end
            arguments (Output)
                stepSize (1,1) double
            end
            stepDirection = -gradient;
            stepSize = obj.m_stepSize;
            if obj.m_maxIterations > 0
                prevCost = obj.evalCost(obj.m_filterSpectrum);
                newCost = obj.evalCost(obj.m_filterSpectrum+stepSize*stepDirection);

                tmp = 0;
                for sIdx = 0:obj.m_numberOfLoudspeakers-1
                    tmp = tmp+real(conj(gradient(:,sIdx+1)).' * stepDirection(:,sIdx+1));
                end
                tmp = tmp*obj.m_filterLength/obj.m_fftSize;
                if newCost > prevCost + obj.m_armijoConstant*stepSize*tmp
                    % Decreasing line search
                    stepSize = obj.decreasingLineSearch(newCost, prevCost, tmp, stepSize, stepDirection);
                else
                    % Increasing line search
                    stepSize = obj.increasingLineSearch(newCost, prevCost,tmp, stepSize, stepDirection);
                end
                if stepSize > obj.m_maxStepSize
                    stepSize = obj.m_maxStepSize;
                elseif stepSize < obj.m_minStepSize
                    stepSize = obj.m_minStepSize;
                end
                obj.m_stepSize = stepSize;
            end
        end
    end
end