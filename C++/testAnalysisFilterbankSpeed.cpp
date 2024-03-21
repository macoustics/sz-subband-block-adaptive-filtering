#include <array>
#include <chrono>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "TfilterbankFast.hpp"

using dataType = double;


const size_t decimationFactor {MY_DECIMATION_FACTOR};
const size_t numberOfSubbands {MY_NUMBER_OF_SUBBANDS};
const size_t prototypeFilterLength {MY_PROTOTYPE_FILTER_LENGTH};
const size_t hopSize {MY_HOP_SIZE};
const size_t numberOfBuffers {MY_NUMBER_OF_BUFFERS};


using inBufferType = std::array<dataType, decimationFactor*hopSize>;
using outBufferType = std::array<std::complex<dataType>, hopSize>;
using outBufferArrayType = std::array<outBufferType, numberOfSubbands/2>;
using prototypeFilterType = std::array<dataType, prototypeFilterLength>;
using analysisFilterbankType = filterbank::TanalysisFilterbank<prototypeFilterType, prototypeFilterLength, inBufferType, outBufferArrayType, dataType, numberOfSubbands, decimationFactor, hopSize>;
using matlabArrayType = matlab::data::TypedArray<dataType>;
using matlabComplexArrayType = matlab::data::TypedArray<std::complex<dataType>>;

using TimeDivision = std::chrono::nanoseconds;
class MexFunction : public matlab::mex::Function {
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;

    std::ostringstream stream;

public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        const matlabArrayType prototypeFilterCoefficients = inputs[0];
        const matlabArrayType inputSamples = inputs[1];
        prototypeFilterType prototypeFilter {};
        inBufferType inputBuffer {};
        outBufferArrayType outputBufferArray {};
        matlabComplexArrayType outputSamples = factory.createArray<std::complex<dataType>>({numberOfSubbands/2,numberOfBuffers*hopSize});
        matlabArrayType outputDurationInNanoseconds = factory.createArray<dataType>({1,1});
        loadFromMatlabArray<matlabArrayType, prototypeFilterType>(prototypeFilterCoefficients, prototypeFilter);

        analysisFilterbankType analysisFilterbank(prototypeFilter);

        auto startTime = std::chrono::system_clock::now();
        for(size_t bIdx = 0; bIdx<numberOfBuffers; bIdx++){
            // Load from inputSamples to inputBuffer
            for(size_t dIdx = 0; dIdx<decimationFactor*hopSize; dIdx++){
                inputBuffer[dIdx] = inputSamples[bIdx*decimationFactor*hopSize + dIdx];
            }
            analysisFilterbank.processBuffer(inputBuffer, outputBufferArray);
            // Transfer outputBuffer to outputSamples
            for(size_t hIdx = 0; hIdx<hopSize; hIdx++){
                for(size_t kIdx = 0; kIdx<numberOfSubbands/2; kIdx++){
                    outputSamples[kIdx][bIdx*hopSize + hIdx] = (outputBufferArray[kIdx])[hIdx];
                }
            }
        }
        auto elapsedTime = std::chrono::duration_cast<TimeDivision>(std::chrono::system_clock::now() - startTime);
        outputDurationInNanoseconds[0] = elapsedTime.count();

        outputs[0] = outputSamples;
        outputs[1] = outputDurationInNanoseconds;
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