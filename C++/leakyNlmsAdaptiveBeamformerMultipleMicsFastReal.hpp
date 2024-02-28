#pragma once
#include <array>
#include <cmath>
#include <boost/circular_buffer.hpp>
#include "gfftRadix4_splitRadix.hpp"

namespace adaptiveFilter{

template<typename dataType, size_t fftSize, size_t numberOfLoudspeakers>
struct settings{
    
};

template<typename realDataType, typename timeDataType, typename frequencyDataType, size_t fftLength, size_t numberOfLoudspeakers, size_t numberOfMicrophones>
struct beamformerSettingsType{
    using sequenceType = std::array<timeDataType, fftLength>;
    using spectrumType = std::array<frequencyDataType, fftLength>;
    using spectraType = std::array<spectrumType, numberOfLoudspeakers>;
    using rirType = std::array<spectraType, numberOfMicrophones>;
    rirType brightRir;
    rirType darkRir;
    spectraType targetFilterSpectra;
    realDataType regularizationFactor;
    realDataType signalEpsilon;
    realDataType minStepSize;
    realDataType maxStepSize;
};

template<typename inBufferType, typename outBufferType, typename realDataType, typename timeDataType, typename frequencyDataType, size_t fftLength, size_t filterLength, size_t numberOfLoudspeakers, size_t numberOfMicrophones>
class beamformerCalculations{
public:
    using sequenceType = std::array<timeDataType, fftLength>;
    using spectrumType = std::array<frequencyDataType, fftLength>;
    using magnitudeSpectrumType = std::array<realDataType, fftLength>;
    using spectraType = std::array<spectrumType, numberOfLoudspeakers>;
    using rirType = std::array<spectraType, numberOfMicrophones>;
    beamformerCalculations(const beamformerSettingsType<realDataType, timeDataType, frequencyDataType, fftLength, numberOfLoudspeakers, numberOfMicrophones>& beamformerSettings)
    : m_inputSpectrum()
    , m_powerSpectrum()
    , m_brightRir(beamformerSettings.brightRir)
    , m_darkRir(beamformerSettings.darkRir)
    , m_filterSpectra()
    , m_targetFilterSpectra(beamformerSettings.targetFilterSpectra)
    , m_signalEpsilon(beamformerSettings.signalEpsilon)
    , m_regularizationFactor(beamformerSettings.regularizationFactor)
    , m_minStepSize(beamformerSettings.minStepSize)
    , m_maxStepSize(beamformerSettings.maxStepSize)
    , m_stepSize(beamformerSettings.maxStepSize/100.0)
    , m_inputBuffer(fftLength)
    {
        realDataType tmpZeroReal {0.0};
        timeDataType tmpZeroTime {0.0};
        for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
            m_inputBuffer.push_back(tmpZeroTime);
            m_powerSpectrum[nIdx] = tmpZeroReal;
        }
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            setSpectrumToZero(m_filterSpectra[sIdx]);
        }
     };
     ~beamformerCalculations() = default;

    boost::circular_buffer<timeDataType> m_inputBuffer;
    spectrumType m_inputSpectrum;
    magnitudeSpectrumType m_powerSpectrum;
    rirType m_brightRir;
    rirType m_darkRir;
    spectraType m_filterSpectra;
    spectraType m_gradient;
    spectraType m_targetFilterSpectra;
    realDataType m_regularizationFactor;
    realDataType m_signalEpsilon;
    realDataType m_gammaSpectrum {0.99};
    realDataType m_stepSize;
    realDataType m_minStepSize;
    realDataType m_maxStepSize;
    realDataType m_armijoConstant {0.0001};
    realDataType m_stepScale {0.1};
    const size_t m_hopSize {fftLength/2};
    const size_t m_maxIterations {2};
    fft::gfft<fftLength, spectrumType, realDataType, true> inverseFft {};
    fft::gfft<fftLength, spectrumType, realDataType, false> fft {};

    void setSpectrumToZero(spectrumType& spectrum){
        frequencyDataType tmpValue {0.0, 0.0};
        for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
            spectrum[nIdx] = tmpValue;
        }
    }

    void mirrorSpectrum(spectrumType& spectrum){
        for(size_t nIdx = 1; nIdx<m_hopSize; nIdx++){
            spectrum[fftLength-nIdx] = std::conj(spectrum[nIdx]);
        }
    }

    void updateInputSpectrum(const inBufferType& inputBuffer){
        for(size_t idx = 0; idx<inputBuffer.size(); idx++){
            m_inputBuffer.push_back(inputBuffer[idx]);
        }
        for(size_t idx = 0; idx<fftLength; idx++){
            m_inputSpectrum[idx].real(m_inputBuffer[idx]);
            m_inputSpectrum[idx].imag(0.0);
        }
        fft.transform(m_inputSpectrum);
        for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
            m_powerSpectrum[nIdx] = m_gammaSpectrum*(m_powerSpectrum[nIdx]) + (1-m_gammaSpectrum)*std::norm(m_inputSpectrum[nIdx]);
        }
    }

    void addGradientStep(const spectraType& inputFilterSpectra, spectraType& outputFilterSpectra){
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                (outputFilterSpectra[sIdx])[nIdx] = (inputFilterSpectra[sIdx])[nIdx] - m_stepSize * ((m_gradient[sIdx])[nIdx]);
            }
        }
    }

    void updateFilters(){
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                (m_filterSpectra[sIdx])[nIdx] -=  m_stepSize * ((m_gradient[sIdx])[nIdx]);
            }
        }
    }

    void applyFilters(outBufferType& outputBuffer){
        spectrumType tmpSpectrum {};
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                tmpSpectrum[nIdx] = m_inputSpectrum[nIdx] * ((m_filterSpectra[sIdx])[nIdx]);
            }
            mirrorSpectrum(tmpSpectrum);
            inverseFft.transform(tmpSpectrum);
            for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
                (outputBuffer[sIdx])[nIdx] = tmpSpectrum[nIdx+m_hopSize].real()/((realDataType)fftLength);
            }
        }
    }

    void calculateBrightError(const spectraType& filterSpectra, const spectraType& rir, spectrumType& errorSpectrum){
        setSpectrumToZero(errorSpectrum);
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                errorSpectrum[nIdx] += ((rir[sIdx])[nIdx]) * (((m_targetFilterSpectra[sIdx])[nIdx]) - ((filterSpectra[sIdx])[nIdx]));
            }
        }
        for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
            errorSpectrum[nIdx] *= m_inputSpectrum[nIdx];
        }
        mirrorSpectrum(errorSpectrum);
        inverseFft.transform(errorSpectrum);
        frequencyDataType tmpZero {0.0, 0.0};
        for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
            errorSpectrum[nIdx] = tmpZero;
        }
        for(size_t nIdx = m_hopSize; nIdx<fftLength; nIdx++){
            errorSpectrum[nIdx] *= 0.5/((realDataType)fftLength);
        }
        fft.transform(errorSpectrum);
    }

    void calculateDarkError(const spectraType& filterSpectra, const spectraType& rir, spectrumType& errorSpectrum){
        setSpectrumToZero(errorSpectrum);
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                errorSpectrum[nIdx] -= ((rir[sIdx])[nIdx]) *  ((filterSpectra[sIdx])[nIdx]);
            }
        }
        for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
            errorSpectrum[nIdx] *= m_inputSpectrum[nIdx];
        }
        mirrorSpectrum(errorSpectrum);
        inverseFft.transform(errorSpectrum);
        frequencyDataType tmpValue {0.0, 0.0};
        for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
            errorSpectrum[nIdx] = tmpValue;
        }
        for(size_t nIdx = m_hopSize; nIdx<fftLength; nIdx++){
            errorSpectrum[nIdx] *= 0.5/((realDataType)fftLength);
        }
        fft.transform(errorSpectrum);
    }

    void updateGradientWithBrightError(const spectraType& filterSpectra, spectraType& gradient){
        spectrumType tmpError {};
        for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
            calculateBrightError(filterSpectra, m_brightRir[mIdx], tmpError);
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                tmpError[nIdx] *= std::conj(m_inputSpectrum[nIdx]);
            }
            for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                    (gradient[sIdx])[nIdx] -= std::conj(((m_brightRir[mIdx])[sIdx])[nIdx]) * (tmpError[nIdx]);
                }
            }
        }
    }

    void updateGradientWithDarkError(const spectraType& filterSpectra, spectraType& gradient){
        spectrumType tmpError {};
        for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
            calculateDarkError(filterSpectra, m_darkRir[mIdx], tmpError);
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                tmpError[nIdx] *= std::conj(m_inputSpectrum[nIdx]);
            }
            for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                    (gradient[sIdx])[nIdx] -= std::conj(((m_darkRir[mIdx])[sIdx])[nIdx]) * (tmpError[nIdx]);
                }
            }
        }
    }

    void normalizeAndRegularize(spectraType& gradient){
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<=m_hopSize; nIdx++){
                (gradient[sIdx])[nIdx] = ((gradient[sIdx])[nIdx]) / (m_powerSpectrum[nIdx] + m_signalEpsilon) + m_regularizationFactor * ((m_filterSpectra[sIdx])[nIdx]);
            }
        }
    }

    void forceLinearConvolution(spectraType& gradient){
        frequencyDataType tmpZero {0.0, 0.0};
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            mirrorSpectrum(gradient[sIdx]);
            inverseFft.transform(gradient[sIdx]);
            for(size_t nIdx = 0; nIdx<filterLength; nIdx++){
                (gradient[sIdx])[nIdx] = (gradient[sIdx])[nIdx]/((realDataType)filterLength);
            }
            for(size_t nIdx = filterLength; nIdx<fftLength; nIdx++){
                (gradient[sIdx])[nIdx] = tmpZero;
            }
            fft.transform(gradient[sIdx]);
        }
    }

    void calculateGradient(const spectraType& filterSpectra, spectraType& gradient){
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            setSpectrumToZero(gradient[sIdx]);
        }
        updateGradientWithBrightError(filterSpectra, gradient);
        updateGradientWithDarkError(filterSpectra, gradient);
        normalizeAndRegularize(gradient);
        forceLinearConvolution(gradient);
    }

    void evaluateCost(const spectraType& filterSpectra, realDataType& cost){
        cost = 0.0;
        spectrumType tmpVector {};
        realDataType tmpValue {};
        for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
            calculateBrightError(filterSpectra, m_brightRir[mIdx], tmpVector);
            for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
                cost += std::norm(tmpVector[nIdx]);
            }
        }
        for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
            calculateDarkError(filterSpectra, m_darkRir[mIdx], tmpVector);
            for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
                cost += std::norm(tmpVector[nIdx]);
            }
        }
        
        tmpValue = 0.0;
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
                tmpValue += std::norm((filterSpectra[sIdx])[nIdx]);
            }
        }
        cost += (m_regularizationFactor/2.0) *tmpValue;
        cost *= 2.0;
    }

    void calculateStepSize(){
        realDataType previousCost {};
        realDataType newCost {};
        spectraType tmpFilterSpectra {};
        evaluateCost(m_filterSpectra, previousCost);
        addGradientStep(m_filterSpectra, tmpFilterSpectra);
        evaluateCost(tmpFilterSpectra, newCost);

        realDataType tmpValue {0.0};
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<m_hopSize; nIdx++){
                tmpValue -= (std::conj((m_gradient[sIdx])[nIdx]) * ((m_gradient[sIdx])[nIdx])).real();
            }
        }

        tmpValue *= 2.0*(realDataType)filterLength/ (realDataType)fftLength;

        if (newCost > previousCost + m_armijoConstant*m_stepSize*tmpValue){
            // Decreasing line search
            size_t lineSearchIdx {0};
            while( (newCost > previousCost + m_armijoConstant*m_stepSize*tmpValue) && (lineSearchIdx < m_maxIterations)){
                lineSearchIdx++;
                m_stepSize *= m_stepScale;
                addGradientStep(m_filterSpectra, tmpFilterSpectra);
                evaluateCost(tmpFilterSpectra, newCost);
            }
        } else{
            // Increasing line search
            size_t lineSearchIdx {0};
            while( (newCost <= previousCost + m_armijoConstant*m_stepSize*tmpValue) && (lineSearchIdx < m_maxIterations)){
                lineSearchIdx++;
                m_stepSize /= m_stepScale;
                addGradientStep(m_filterSpectra, tmpFilterSpectra);
                evaluateCost(tmpFilterSpectra, newCost);
            }
            // Choose the last stepSize satisfying the Armijo condition
            if(lineSearchIdx > 0){
                m_stepSize *= m_stepScale;
            }
        }
        if(m_stepSize > m_maxStepSize){
            m_stepSize = m_maxStepSize;
        }else if(m_stepSize < m_minStepSize){
            m_stepSize = m_minStepSize;
        }
    }
};

template<typename inBufferType, typename outBufferType, typename realDataType, typename timeDataType, typename frequencyDataType, size_t fftLength, size_t filterLength, size_t numberOfLoudspeakers, size_t numberOfMicrophones>
class beamformer{
public:
    beamformer(const beamformerSettingsType<realDataType, timeDataType, frequencyDataType, fftLength, numberOfLoudspeakers, numberOfMicrophones>& beamformerSettings)
    : m_beamformerCalculations(beamformerSettings)
    {};
    ~beamformer() = default;

    void processInputBuffer(const inBufferType& inputBuffer, outBufferType& outputBuffer){
        m_beamformerCalculations.updateInputSpectrum(inputBuffer);
        m_beamformerCalculations.calculateGradient(m_beamformerCalculations.m_filterSpectra, m_beamformerCalculations.m_gradient);
        m_beamformerCalculations.calculateStepSize();
        m_beamformerCalculations.updateFilters();
        m_beamformerCalculations.applyFilters(outputBuffer);
    }
private:
    using beamformerCalculationsType = beamformerCalculations<inBufferType, outBufferType, realDataType, timeDataType, frequencyDataType, fftLength, filterLength, numberOfLoudspeakers, numberOfMicrophones>;
    beamformerCalculationsType m_beamformerCalculations;
};



} // end namespace