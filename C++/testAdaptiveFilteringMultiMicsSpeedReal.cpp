#include <array>
#include <chrono>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "leakyNlmsAdaptiveBeamformerMultipleMicsFastReal.hpp"

using dataType = double;
using cDataType = std::complex<dataType>;
const size_t filterLength {MY_FILTER_LENGTH};
const size_t numberOfLoudspeakers {MY_NUMBER_OF_LOUDSPEAKERS};
const size_t numberOfMicrophones {MY_NUMBER_OF_MICROPHONES};
const size_t numberOfBuffers {MY_NUMBER_OF_BUFFERS}; // 
const size_t fftLength {MY_FFT_LENGTH};
const size_t hopSize {MY_HOP_SIZE};


using inBufferType = std::array<dataType, hopSize>;
using outBufferType = std::array<inBufferType, numberOfLoudspeakers>;
using beamformerSettingsType = adaptiveFilter::beamformerSettingsType<dataType, dataType, cDataType, fftLength, numberOfLoudspeakers, numberOfMicrophones>;
using beamformerType = adaptiveFilter::beamformer<inBufferType, outBufferType, dataType, dataType, cDataType, fftLength, filterLength, numberOfLoudspeakers, numberOfMicrophones>;

using matlabArrayType = matlab::data::TypedArray<dataType>;
using matlabComplexArrayType = matlab::data::TypedArray<cDataType>;
using matlabBoolArrayType = matlab::data::TypedArray<bool>;

using TimeDivision = std::chrono::nanoseconds;

class MexFunction : public matlab::mex::Function {
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;

    std::ostringstream stream;

public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        const matlabArrayType inputSamples = (inputs[0]);
        const matlabComplexArrayType matlabBrightRir = (inputs[1]);
        const matlabComplexArrayType matlabDarkRir = (inputs[2]);
        const matlabComplexArrayType matlabTargetFilterSpectra = (inputs[3]);
        const matlabArrayType matlabSignalEpsilon = inputs[4];
        const matlabArrayType matlabRegParam = inputs[5];
        const matlabArrayType matlabStepRange = inputs[6];
        const size_t numberOfRepeats {1};

        matlabArrayType outputSamples = factory.createArray<dataType>({hopSize*numberOfBuffers,numberOfLoudspeakers});
        matlabArrayType outputDurationInNanoseconds = factory.createArray<dataType>({1,1});
        
        inBufferType tmpInputBuffer {};

        outBufferType tmpOutputBuffer {};
        beamformerSettingsType beamformerSettings {};
        

        for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
            for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
                    ((beamformerSettings.brightRir[mIdx])[sIdx])[nIdx] = (cDataType)(matlabBrightRir[nIdx][sIdx][mIdx]);
                    ((beamformerSettings.darkRir[mIdx])[sIdx])[nIdx] = (cDataType)(matlabDarkRir[nIdx][sIdx][mIdx]);
                }
            }
        }
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
                (beamformerSettings.targetFilterSpectra[sIdx])[nIdx] = (cDataType)(matlabTargetFilterSpectra[nIdx][sIdx]);
            }
        }
        beamformerSettings.regularizationFactor = matlabRegParam[0];
        beamformerSettings.signalEpsilon = matlabSignalEpsilon[0];
        beamformerSettings.minStepSize = matlabStepRange[0];
        beamformerSettings.maxStepSize = matlabStepRange[1];
        beamformerType beamformer(beamformerSettings);

        auto startTime = std::chrono::system_clock::now();
        for(size_t repeatIdx = 0; repeatIdx<numberOfRepeats; repeatIdx++){
            for(size_t bIdx = 0; bIdx<numberOfBuffers; bIdx++){
                // Each block is of size decimationFactor * hopSize
                for(size_t hIdx = 0; hIdx<hopSize; hIdx++){
                    // Apply the analysisFilterbank hopSize times
                    tmpInputBuffer[hIdx] = inputSamples[bIdx*hopSize + hIdx];
                }
                beamformer.processInputBuffer(tmpInputBuffer,tmpOutputBuffer);
                
                // Synthesize each of the loudspeaker signals
                for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                    for(size_t idx = 0; idx<hopSize; idx++){
                        outputSamples[bIdx*hopSize + idx][sIdx] = (tmpOutputBuffer[sIdx])[idx];
                    }
                }
            }
        }
        auto elapsedTime = std::chrono::duration_cast<TimeDivision>(std::chrono::system_clock::now() - startTime);
        outputDurationInNanoseconds[0] = elapsedTime.count()/numberOfRepeats;

        outputs[0] = std::move(outputSamples);
        outputs[1] = std::move(outputDurationInNanoseconds);
        
        
    }      
    template<typename matlabType, typename cppType>
    void loadFromMatlabArray(const matlabType& inputArray, cppType& outputArray){
        for(size_t idx=0; idx<inputArray.getNumberOfElements(); idx++){
            outputArray[idx] = inputArray[idx];
        }
    }

    template<typename matlabType, typename cppType>
    void saveToMatlabArray(const cppType& inputArray, matlabType& outputArray){
        for(size_t idx=0; idx<inputArray.size(); idx++){
            outputArray[idx] = inputArray[idx];
        }
    }

    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};