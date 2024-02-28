clear all; close all; clc
% 
buildTask = 'analysisFilterbank';
switch buildTask
    case 'analysisFilterbank'
        eval(['mex -f mexopts.xml testAnalysisFilterbankSpeed.cpp -IC:\boost_1_82_0\'])
    case 'real'
        eval(['mex -f mexopts.xml testAdaptiveFilteringMultiMicsSpeedReal.cpp -IC:\boost_1_82_0\'])
    case 'mics'
        eval(['mex -f mexopts.xml testSubbandAdaptiveFilteringMultipleMics.cpp -IC:\boost_1_82_0\'])
    case 'parallel'
        eval(['mex -f mexopts.xml testSubbandAdaptiveFilteringMicsParallel.cpp -IC:\boost_1_82_0\'])
end