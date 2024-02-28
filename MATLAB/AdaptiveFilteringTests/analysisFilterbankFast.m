classdef analysisFilterbankFast < handle
    properties
        m_decimationFactor = 0;
        m_numberOfSubbands = 0;
        m_prototypeFilterLength = 0;
        m_prototypeFilterLengthModSubbands = 0;
        m_inputBuffer = [];
        m_numberOfChannelsMultiplier = 0;
        m_Tmatrix = [];
        m_prototypeFilter = [];
        m_D1 = [];
        m_D2 = [];
    end
    methods 
        function analysisFilterbankFast = analysisFilterbankFast(prototypeFilter, numberOfSubbands, decimationFactor)
            %% Initialize polyphase matrix
            prototypeFilterLength = length(prototypeFilter);
            for idx = 0:prototypeFilterLength-1
                if mod(floor(idx/numberOfSubbands),2)
                    prototypeFilter(idx+1) = -prototypeFilter(idx+1);
                end
            end
            numberOfChannelsMultiplier = ceil(prototypeFilterLength/numberOfSubbands);
            analysisFilterbankFast.m_numberOfChannelsMultiplier = numberOfChannelsMultiplier;
            analysisFilterbankFast.m_decimationFactor = decimationFactor;
            analysisFilterbankFast.m_numberOfSubbands = numberOfSubbands;
            analysisFilterbankFast.m_prototypeFilter = prototypeFilter;
            analysisFilterbankFast.m_prototypeFilterLength = prototypeFilterLength;

            %% Create the input-delay buffer
            analysisFilterbankFast.m_inputBuffer = zeros(prototypeFilterLength,1);
            timeOffset = -(prototypeFilterLength-1)/2;
            channelOffset = 0.5;
            time = (0:numberOfSubbands-1)';
            channels = (0:numberOfSubbands/2-1)';
            
            analysisFilterbankFast.m_D1 = exp(1j*2*pi*(channels*timeOffset/numberOfSubbands));
            analysisFilterbankFast.m_D2 = exp(1j*2*pi*(channelOffset*(time+timeOffset)/numberOfSubbands));
%             T = exp(1j*2*pi/numberOfSubbands*(channels + channel_offset) * (time' + time_offset));
%             analysisFilterbankFast.m_Tmatrix = T;
            analysisFilterbankFast.m_prototypeFilterLengthModSubbands = mod(prototypeFilterLength, numberOfSubbands);
        end

        function [obj] = updateInputDelayBuffer(obj, inputBuffer)
            % Note that it is expected that the input buffer is arranged
            % from oldest sample to newest sample
            arguments (Input)
                obj (1,1) analysisFilterbankFast
                inputBuffer (:,1) double
            end
            arguments (Output)
                obj (1,1) analysisFilterbankFast
            end
            for dIdx = 0:obj.m_decimationFactor-1
                obj.m_inputBuffer = [inputBuffer(dIdx+1); obj.m_inputBuffer(1:end-1)];
            end
        end

        function [filteredSamples, obj] = applyPolyphaseFilter(obj)
            arguments (Input)
                obj (1,1) analysisFilterbankFast
            end
            arguments (Output)
                filteredSamples (:,1) double
                obj (1,1) analysisFilterbankFast
            end
            filteredSamples = obj.m_prototypeFilter .* obj.m_inputBuffer;
        end

        function [outputSamples] = determineOutputSamples(obj, filteredSamples)
            arguments (Input)
                obj (1,1) analysisFilterbankFast
                filteredSamples (:,1) double
            end
            arguments (Output)
                outputSamples (:,1) double
            end
            outputSamples = zeros(obj.m_numberOfSubbands,1);
            for jIdx = 0:obj.m_numberOfChannelsMultiplier-2
                idx = jIdx*obj.m_numberOfSubbands + (1:obj.m_numberOfSubbands);
                outputSamples = outputSamples + filteredSamples(idx);
            end
            for idx = 0:obj.m_prototypeFilterLengthModSubbands-1
                outputSamples(idx+1) = outputSamples(idx+1) + filteredSamples((obj.m_numberOfChannelsMultiplier-1)*obj.m_numberOfSubbands + idx+1);
            end
%             outputSamples = obj.m_Tmatrix*outputSamples;
            tmpSamples = obj.m_numberOfSubbands *ifft(obj.m_D2.*outputSamples);
            outputSamples =  obj.m_D1.*tmpSamples(1:obj.m_numberOfSubbands/2);
        end

        function [outputBuffer, obj] = processInputBuffer(obj, inputBuffer)
            arguments (Input)
                obj (1,1) analysisFilterbankFast
                inputBuffer (:,1) double
            end
            arguments (Output)
                outputBuffer (:,1) double
                obj (1,1) analysisFilterbankFast
            end
            % Update the input buffer
            obj.updateInputDelayBuffer(inputBuffer);
            % Calculate the output samples
            outputBuffer = obj.determineOutputSamples(obj.applyPolyphaseFilter);
        end
    end
end