#include <array>
#include <chrono>
#include <omp.h>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "leakyNlmsAdaptiveBeamformerMultipleMicsFast.hpp"
#include "TfilterbankFast.hpp"

using dataType = double;
using cDataType = std::complex<dataType>;
const size_t decimationFactor {4};
const size_t numberOfSubbands {8};
const size_t prototypeFilterLength {32};
const size_t filterLength {75};
const size_t numberOfLoudspeakers {7};
const size_t numberOfMicrophones {3};
const size_t numberOfBuffers {1406}; 
const size_t fftLength {512};
const size_t hopSize {256};


using inBufferType = std::array<dataType, decimationFactor*hopSize>;
using adaptiveFilteringInBufferType = std::array<cDataType, hopSize>;
using analysisBufferArrayType = std::array<adaptiveFilteringInBufferType, numberOfSubbands/2>;
using adaptiveFilteringOutBufferType = std::array<adaptiveFilteringInBufferType, numberOfLoudspeakers>;
using parallelAnalysisBufferArrayType = std::array<analysisBufferArrayType, numberOfLoudspeakers>;
using synthesisBufferArrayType = std::array<adaptiveFilteringOutBufferType, numberOfSubbands/2>;

using outBufferType = std::array<inBufferType, numberOfLoudspeakers>;
using beamformerSettingsType = adaptiveFilter::beamformerSettingsType<dataType, cDataType, cDataType, fftLength, numberOfLoudspeakers, numberOfMicrophones>;
using beamformerSettingsArrayType = std::array<beamformerSettingsType, numberOfSubbands/2>;
using beamformerType = adaptiveFilter::beamformer<adaptiveFilteringInBufferType, adaptiveFilteringOutBufferType, dataType, cDataType, cDataType, fftLength, filterLength, numberOfLoudspeakers, numberOfMicrophones>;
using beamformerPtrType = std::shared_ptr<beamformerType>;
using beamformerPtrArrayType = std::array<beamformerPtrType, numberOfSubbands/2>;

using prototypeFilterType = std::array<dataType, prototypeFilterLength>;
using analysisFilterbankType = filterbank::TanalysisFilterbank<prototypeFilterType, prototypeFilterLength, inBufferType, analysisBufferArrayType, dataType, numberOfSubbands, decimationFactor, hopSize>;
using synthesisFilterbankType = filterbank::TsynthesisFilterbank<prototypeFilterType, prototypeFilterLength, analysisBufferArrayType, inBufferType, dataType, numberOfSubbands, decimationFactor, hopSize>;
using synthesisFilterbankPtrType = std::shared_ptr<synthesisFilterbankType>;
using synthesisFilterbankPtrArrayType = std::array<synthesisFilterbankPtrType, numberOfLoudspeakers>;

using booleanGainArray = std::array<bool, numberOfSubbands/2>;

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
        const matlabArrayType prototypeFilterCoefficients = std::move(inputs[0]);
        const matlabArrayType inputSamples = std::move(inputs[1]);
        const matlabComplexArrayType matlabBrightRir = std::move(inputs[2]);
        const matlabComplexArrayType matlabDarkRir = std::move(inputs[3]);
        const matlabComplexArrayType matlabTargetFilterSpectra = std::move(inputs[4]);
        const matlabArrayType matlabSignalEpsilon = inputs[5];
        const matlabArrayType matlabRegParam = inputs[6];
        const matlabArrayType matlabStepRange = inputs[7];
        const matlabBoolArrayType matlabBinaryGains = inputs[8];

        matlabArrayType outputSamples = factory.createArray<dataType>({hopSize*decimationFactor*numberOfBuffers,numberOfLoudspeakers});
        matlabArrayType outputDurationInNanoseconds = factory.createArray<dataType>({1,1});
        
        inBufferType tmpInputBuffer {};
        outBufferType tmpOutputBuffer {};
        prototypeFilterType prototypeFilter {};

        analysisBufferArrayType analysisBufferArray {};
        parallelAnalysisBufferArrayType tmpSynthesisBuffer {};
        synthesisBufferArrayType synthesisBufferArray {};
        cDataType tmpZero {0.0, 0.0};
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            for(size_t cIdx = 0; cIdx<numberOfSubbands/2; cIdx++){
                for(size_t nIdx = 0; nIdx<hopSize; nIdx++){
                    ((tmpSynthesisBuffer[sIdx])[cIdx])[nIdx] = tmpZero;
                }
            }
        }

        beamformerSettingsArrayType beamformerSettingsArray {};
        beamformerPtrArrayType beamformerPtrArray;
        booleanGainArray binaryGainsArray {};

        for(size_t cIdx = 0; cIdx<numberOfSubbands/2; cIdx++){
            for(size_t mIdx = 0; mIdx<numberOfMicrophones; mIdx++){
                for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                    for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
                        ((beamformerSettingsArray[cIdx].brightRir[mIdx])[sIdx])[nIdx] = (cDataType)(matlabBrightRir[nIdx][sIdx][mIdx][cIdx]);
                        ((beamformerSettingsArray[cIdx].darkRir[mIdx])[sIdx])[nIdx] = (cDataType)(matlabDarkRir[nIdx][sIdx][mIdx][cIdx]);
                    }
                }
            }
            for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                for(size_t nIdx = 0; nIdx<fftLength; nIdx++){
                    (beamformerSettingsArray[cIdx].targetFilterSpectra[sIdx])[nIdx] = (cDataType)(matlabTargetFilterSpectra[nIdx][sIdx][cIdx]);
                }
            }

            beamformerSettingsArray[cIdx].regularizationFactor = matlabRegParam[cIdx];
            beamformerSettingsArray[cIdx].signalEpsilon = matlabSignalEpsilon[cIdx];
            beamformerSettingsArray[cIdx].minStepSize = matlabStepRange[0][cIdx];
            beamformerSettingsArray[cIdx].maxStepSize = matlabStepRange[1][cIdx];
            beamformerPtrArray[cIdx] = std::make_shared<beamformerType>(beamformerSettingsArray[cIdx]);
        }

        loadFromMatlabArray<matlabArrayType, prototypeFilterType>(prototypeFilterCoefficients, prototypeFilter);
        loadFromMatlabArray<matlabBoolArrayType, booleanGainArray>(matlabBinaryGains, binaryGainsArray);
        analysisFilterbankType analysisFilterbank(prototypeFilter);
        synthesisFilterbankPtrArrayType synthesisFilterbankPtrArray;
        for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
            synthesisFilterbankPtrArray[sIdx] = std::make_shared<synthesisFilterbankType>(prototypeFilter);
        }

        auto startTime = std::chrono::system_clock::now();
        for(size_t bIdx = 0; bIdx<numberOfBuffers; bIdx++){
            // Each block is of size decimationFactor * hopSize
            for(size_t dIdx = 0; dIdx<decimationFactor*hopSize; dIdx++){
                // Apply the analysisFilterbank hopSize times
                tmpInputBuffer[dIdx] = inputSamples[bIdx*hopSize*decimationFactor + dIdx];
            }
            analysisFilterbank.processBuffer(tmpInputBuffer, analysisBufferArray);
            // Process the adaptive filters in each subband
            #pragma omp parallel for
            for(size_t cIdx = 0; cIdx<numberOfSubbands/2; cIdx++){
                //#pragma omp single
                {
                    if(binaryGainsArray[cIdx]){
                        beamformerPtrArray[cIdx]->processInputBuffer(analysisBufferArray[cIdx],synthesisBufferArray[cIdx]);
                    }
                }
            }
            // Synthesize each of the loudspeaker signals
            #pragma omp parallel for
            for(size_t sIdx = 0; sIdx<numberOfLoudspeakers; sIdx++){
                //#pragma omp single
                {
                    for(size_t cIdx = 0; cIdx<numberOfSubbands/2; cIdx++){
                        if(binaryGainsArray[cIdx]){
                            for(size_t nIdx = 0; nIdx<hopSize; nIdx++){
                                ((tmpSynthesisBuffer[sIdx])[cIdx])[nIdx] = ((synthesisBufferArray[cIdx])[sIdx])[nIdx];
                            }
                        }
                    }
                    synthesisFilterbankPtrArray[sIdx]->processBuffer(tmpSynthesisBuffer[sIdx], tmpOutputBuffer[sIdx]);
                    for(size_t idx = 0; idx<decimationFactor*hopSize; idx++){
                        outputSamples[bIdx*hopSize*decimationFactor + idx][sIdx] = (tmpOutputBuffer[sIdx])[idx];
                    }
                }
            }
        }
        auto elapsedTime = std::chrono::duration_cast<TimeDivision>(std::chrono::system_clock::now() - startTime);
        outputDurationInNanoseconds[0] = elapsedTime.count();

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