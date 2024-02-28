classdef synthesisFilterbankFast < handle
    properties
        m_decimationFactor = 0;
        m_numberOfSubbands = 0;
        m_prototypeFilterLength = 0;
        m_prototypeFilterLengthModSubbands = 0;
        m_outputBuffer = [];
        m_numberOfChannelsMultiplier = 0;
        m_Tmatrix = [];
        m_prototypeFilter = [];
        m_D1conj = [];
        m_D2conj = [];
    end
    methods 
        function synthesisFilterbankFast = synthesisFilterbankFast(prototypeFilter, numberOfSubbands, decimationFactor)
            %% Initialize polyphase matrix
            prototypeFilterLength = length(prototypeFilter);
            for idx = 0:prototypeFilterLength-1
                if mod(floor(idx/numberOfSubbands),2)
                    prototypeFilter(idx+1) = -prototypeFilter(idx+1);
                end
            end
            numberOfChannelsMultiplier = ceil(prototypeFilterLength/numberOfSubbands);
            synthesisFilterbankFast.m_numberOfChannelsMultiplier = numberOfChannelsMultiplier;
            synthesisFilterbankFast.m_decimationFactor = decimationFactor;
            synthesisFilterbankFast.m_numberOfSubbands = numberOfSubbands;
            synthesisFilterbankFast.m_prototypeFilter = prototypeFilter;
            synthesisFilterbankFast.m_prototypeFilterLength = prototypeFilterLength;

            %% Create the input-delay buffer
            synthesisFilterbankFast.m_outputBuffer = zeros(prototypeFilterLength,1);
            timeOffset = -(prototypeFilterLength-1)/2;
            channelOffset = 0.5;
            time = (0:numberOfSubbands-1)';
            channels = (0:numberOfSubbands/2-1)';
            
            synthesisFilterbankFast.m_D1conj = exp(-1j*2*pi*(channels*timeOffset/numberOfSubbands));
            synthesisFilterbankFast.m_D2conj = exp(-1j*2*pi*(channelOffset*(time+timeOffset)/numberOfSubbands));
%             T = exp(1j*2*pi/numberOfSubbands*(channels + channel_offset) * (time' + time_offset));
%             synthesisFilterbankFast.m_Tmatrix = T;
            synthesisFilterbankFast.m_prototypeFilterLengthModSubbands = mod(prototypeFilterLength, numberOfSubbands);
        end

        function [filteredSamples, obj] = applyPolyphaseFilter(obj, inputSamples)
            arguments (Input)
                obj (1,1) synthesisFilterbankFast
                inputSamples (:,1) double
            end
            arguments (Output)
                filteredSamples (:,1) double
                obj (1,1) synthesisFilterbankFast
            end
            tmpSamples = obj.m_D1conj.*inputSamples;
            tmpSamples = 2*real(obj.m_D2conj.*fft(tmpSamples,obj.m_numberOfSubbands));
            filteredSamples = zeros(obj.m_prototypeFilterLength,1);
            filteredSamples(1:(obj.m_numberOfChannelsMultiplier-1)*obj.m_numberOfSubbands) = repmat(tmpSamples,obj.m_numberOfChannelsMultiplier-1,1);
            for idx = 0:obj.m_prototypeFilterLengthModSubbands
                filteredSamples((obj.m_numberOfChannelsMultiplier-1)*obj.m_numberOfSubbands + idx + 1) = tmpSamples(idx+1);
            end

            filteredSamples = filteredSamples.*obj.m_prototypeFilter;
        end

        function [outputBuffer, obj] = determineOutputSamples(obj, filteredSamples)
            arguments (Input)
                obj (1,1) synthesisFilterbankFast
                filteredSamples (:,1) double
            end
            arguments (Output)
                outputBuffer (:,1) double
                obj (1,1) synthesisFilterbankFast
            end
            % Calculate output samples
            obj.m_outputBuffer = obj.m_outputBuffer + filteredSamples;
            outputBuffer = zeros(obj.m_decimationFactor,1);
            for idx=0:obj.m_decimationFactor-1
                outputBuffer(idx+1) = obj.m_outputBuffer(obj.m_prototypeFilterLength-idx);
            end
            
            % Update the output delay line
            obj.m_outputBuffer(obj.m_decimationFactor+1:end) = obj.m_outputBuffer(1:end-obj.m_decimationFactor);
            obj.m_outputBuffer(1:obj.m_decimationFactor) = zeros(obj.m_decimationFactor,1);
        end

        function [outputBuffer, obj] = processInputBuffer(obj, inputBuffer)
            arguments (Input)
                obj (1,1) synthesisFilterbankFast
                inputBuffer (:,1) double
            end
            arguments (Output)
                outputBuffer (:,1) double
                obj (1,1) synthesisFilterbankFast
            end
            outputBuffer = obj.determineOutputSamples(obj.applyPolyphaseFilter(inputBuffer));
        end
    end
end