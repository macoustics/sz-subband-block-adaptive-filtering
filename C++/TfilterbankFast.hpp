#pragma once
#include <array>
#include <complex>
#include <boost/circular_buffer.hpp>
#include "gfftRadix4_splitRadix.hpp"

namespace filterbank{
template <typename T>
constexpr T intCeilDiv(const T& x, const T& y){
    return (x + y -1)/y;
}


template <typename prototypeFilterType, size_t prototypeFilterLength, typename inBufferType, typename outBufferType, typename dataType, size_t numberOfSubbands, size_t decimationFactor>
struct TfilterbankSettings{
    using complexArrayType = std::array<std::complex<dataType>,numberOfSubbands>;
    using realArrayType = std::array<dataType,numberOfSubbands>;
    using ringBufferType = boost::circular_buffer<dataType>;
    const dataType subbandChannelOffset {0.5};
    const dataType timeOffset {-(dataType)(prototypeFilterLength-1)/2.0};
    const dataType pi {std::acos((dataType)-1)};
    const size_t channelMultiplier {prototypeFilterLength/numberOfSubbands};
    const size_t prototypeFilterLengthModSubbands {prototypeFilterLength % numberOfSubbands};
    std::array<dataType, prototypeFilterLength> prototypeFilter;
    std::array<std::complex<dataType>, numberOfSubbands/2> D1;
    std::array<std::complex<dataType>, numberOfSubbands> D2;
    std::unique_ptr<ringBufferType> ringBufferPtr;
};


template <typename prototypeFilterType, size_t prototypeFilterLength, typename inBufferType, typename outBufferType, typename dataType, size_t numberOfSubbands, size_t decimationFactor, size_t hopSize>
class TanalysisFilterbank {
public:
    using SettingsType = TfilterbankSettings<prototypeFilterType, prototypeFilterLength, inBufferType, outBufferType, dataType, numberOfSubbands, decimationFactor>;
    using cDataType = std::complex<dataType>;
    TanalysisFilterbank(const prototypeFilterType& prototypeFilter) {
        fftPtr = std::make_shared<fftType>();
        m_settings.ringBufferPtr = std::make_unique<typename SettingsType::ringBufferType>(prototypeFilterLength);
        // Populate D1
        for(size_t idx = 0; idx<numberOfSubbands/2; idx++){
            m_settings.D1[idx] = std::exp(cDataType(0, (dataType)idx * 2.0 * m_settings.pi * m_settings.timeOffset / (dataType)numberOfSubbands));
        }
        for(size_t idx = 0; idx<numberOfSubbands; idx++){
            m_settings.D2[idx] = std::exp(cDataType(0, ((dataType)idx + m_settings.timeOffset) * 2.0 * m_settings.pi * m_settings.subbandChannelOffset / (dataType)numberOfSubbands));
        }
        
        for(size_t idx = 0; idx<prototypeFilterLength; idx++){
            if ((idx/numberOfSubbands) % 2){
                m_settings.prototypeFilter[idx] = -prototypeFilter[idx];
            } else{
                m_settings.prototypeFilter[idx] = prototypeFilter[idx];
            }
            m_settings.ringBufferPtr->push_back(0.0);
        }
    }

    void processBuffer(const inBufferType& inputBuffer, outBufferType& outputBuffer){
        typename SettingsType::realArrayType tmpRealBuffer {};
        typename SettingsType::complexArrayType tmpComplexBuffer {};
        for (size_t hIdx = 0; hIdx<hopSize; hIdx++){
            for (size_t dIdx = 0; dIdx<decimationFactor; dIdx++){
                tmpRealBuffer[dIdx] = inputBuffer[hIdx*decimationFactor + dIdx];
            }
            updateInputBuffer(tmpRealBuffer);
            applyPrototypeFilter(tmpRealBuffer);
            for(size_t idx = 0; idx<numberOfSubbands; idx++){
                tmpComplexBuffer[idx] = tmpRealBuffer[idx] * m_settings.D2[idx];
            }
            fftPtr->transform(tmpComplexBuffer);
            for(size_t idx = 0; idx<numberOfSubbands/2; idx++){
                (outputBuffer[idx])[hIdx] = m_settings.D1[idx] * tmpComplexBuffer[idx];
            }
        }
    }

    ~TanalysisFilterbank() = default;
protected:
    using fftType = fft::gfft<numberOfSubbands, typename SettingsType::complexArrayType, dataType, true>;
    std::shared_ptr<fftType> fftPtr;
    SettingsType m_settings {};

    void updateInputBuffer(const typename SettingsType::realArrayType& inputBuffer){
        for(size_t dIdx = 0; dIdx<decimationFactor; dIdx++){
            m_settings.ringBufferPtr->push_front(inputBuffer[dIdx]);
        }
    }

    void applyPrototypeFilter(typename SettingsType::realArrayType& outputBuffer){
        size_t tmpIdx {};
        for(size_t cIdx=0; cIdx<numberOfSubbands; cIdx++){
            outputBuffer[cIdx] = 0.0;
        }
        for(size_t jIdx = 0; jIdx<m_settings.channelMultiplier; jIdx++){
            for(size_t cIdx = 0; cIdx<numberOfSubbands; cIdx++){
                tmpIdx = jIdx*numberOfSubbands + cIdx;
                outputBuffer[cIdx] += (*m_settings.ringBufferPtr)[tmpIdx] * m_settings.prototypeFilter[tmpIdx];
            }
        }
        for(size_t cIdx = 0; cIdx<m_settings.prototypeFilterLengthModSubbands; cIdx++){
            tmpIdx = m_settings.channelMultiplier*numberOfSubbands;
            outputBuffer[cIdx] += (*m_settings.ringBufferPtr)[tmpIdx+cIdx] * m_settings.prototypeFilter[tmpIdx+cIdx];
        }
    }
};


template <typename prototypeFilterType, size_t prototypeFilterLength, typename inBufferType, typename outBufferType, typename dataType, size_t numberOfSubbands, size_t decimationFactor, size_t hopSize>
class TsynthesisFilterbank {
public:
    using SettingsType = TfilterbankSettings<prototypeFilterType, prototypeFilterLength, inBufferType, outBufferType, dataType, numberOfSubbands, decimationFactor>;
    using cDataType = std::complex<dataType>;
    TsynthesisFilterbank(const prototypeFilterType& prototypeFilter){
        fftPtr = std::make_shared<fftType>();
        m_settings.ringBufferPtr = std::make_unique<typename SettingsType::ringBufferType>(prototypeFilterLength);
        // Populate D1
        for(size_t idx = 0; idx<numberOfSubbands/2; idx++){
            m_settings.D1[idx] = std::exp(cDataType(0, -(dataType)idx * 2.0 * m_settings.pi * m_settings.timeOffset / (dataType)numberOfSubbands));
        }
        for(size_t idx = 0; idx<numberOfSubbands; idx++){
            m_settings.D2[idx] = std::exp(cDataType(0, -((dataType)idx + m_settings.timeOffset) * 2.0 * m_settings.pi * m_settings.subbandChannelOffset / (dataType)numberOfSubbands));
        }
        for(size_t idx = 0; idx<prototypeFilterLength; idx++){
            if ((idx/numberOfSubbands) % 2){
                m_settings.prototypeFilter[idx] = -prototypeFilter[idx];
            } else{
                m_settings.prototypeFilter[idx] = prototypeFilter[idx];
            }
            m_settings.ringBufferPtr->push_back(0.0);
        }
    }

    void processBuffer(const inBufferType& inputBuffer, outBufferType& outputBuffer){
        for(size_t hIdx = 0; hIdx<hopSize; hIdx++){
            updateOutputBuffer(inputBuffer, hIdx);
            determineOutputSamples(outputBuffer, hIdx);
        }
    }

    ~TsynthesisFilterbank() = default;
protected:
    using fftType = fft::gfft<numberOfSubbands, typename SettingsType::complexArrayType, dataType, false>;
    std::shared_ptr<fftType> fftPtr;
    SettingsType m_settings {};

    void updateOutputBuffer(const inBufferType& inputBuffer, const size_t& hIdx){
        size_t tmpIdx {};
        cDataType tmpComplexValue {};
        typename SettingsType::realArrayType tmpRealBuffer {};
        typename SettingsType::complexArrayType tmpComplexBuffer {};
        for(size_t idx = 0; idx<numberOfSubbands/2; idx++){
            tmpComplexBuffer[idx] = idx < numberOfSubbands/2 ? 2.0 * m_settings.D1[idx] * (inputBuffer[idx])[hIdx] : cDataType(0,0);
        }
        fftPtr->transform(tmpComplexBuffer);
        for(size_t idx = 0; idx<numberOfSubbands; idx++){
             tmpComplexBuffer[idx] *= m_settings.D2[idx];
        }
        for(size_t jIdx = 0; jIdx<m_settings.channelMultiplier; jIdx++){
            for(size_t cIdx = 0; cIdx<numberOfSubbands; cIdx++){
                tmpIdx = jIdx*numberOfSubbands + cIdx;
                (*m_settings.ringBufferPtr)[tmpIdx] += m_settings.prototypeFilter[tmpIdx] * tmpComplexBuffer[cIdx].real();
            }
        }
        for(size_t cIdx = 0; cIdx<m_settings.prototypeFilterLengthModSubbands; cIdx++){
            tmpIdx = m_settings.channelMultiplier*numberOfSubbands;
            (*m_settings.ringBufferPtr)[tmpIdx+cIdx] += m_settings.prototypeFilter[tmpIdx+cIdx] * tmpComplexBuffer[cIdx].real();
        }
    }

    void determineOutputSamples(outBufferType& outputBuffer, const size_t hIdx){
        for(size_t idx = 0; idx < decimationFactor; idx++){
            outputBuffer[hIdx*decimationFactor + idx] = (*m_settings.ringBufferPtr)[prototypeFilterLength-1-idx];
        }
        for(size_t idx=0; idx<decimationFactor; idx++){
            m_settings.ringBufferPtr->push_front(0.0);
        }
    }

};

} // Namespace end