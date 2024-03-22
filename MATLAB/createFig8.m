clear all; close all; clc

task = 'plotMisalignment';
numberOfChannels = 8;
oversamplingFactor = 5/4;
decimationFactor = round(numberOfChannels/oversamplingFactor);
prototypeFilterLength = 16*numberOfChannels;
setupEnvironment = 'Simulation';
evaluationEnvironment = 'anechoic';
numberOfRepeats = 1;

%% Create Fig. 8 a
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, false, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Fig. 8 (a)');

%% Create Fig. 8 b
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, false, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Fig. 8 (b)');

%% Create Fig. 8 c
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, false, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Fig. 8 (c)');

%% Create Fig. 8 d
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, false, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Fig. 8 (d)');

%% Create Fig. 8 e (REQUIRES MUSIC SIGNAL)
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, true, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Fig. 8 (e)');

%% Create Fig. 8 f (REQUIRES MUSIC SIGNAL)
runAdaptiveProcessing(task, 'matlabFullrate', numberOfChannels, decimationFactor, prototypeFilterLength, true, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Fig. 8 (f)');

%% Create Fig. 8 g (REQUIRES MUSIC SIGNAL)
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, true, setupEnvironment, evaluationEnvironment, numberOfRepeats, 0, true);
title('Fig. 8 (g)');

%% Create Fig. 8 h (REQUIRES MUSIC SIGNAL)
runAdaptiveProcessing(task, 'matlabSubband', numberOfChannels, decimationFactor, prototypeFilterLength, true, setupEnvironment, evaluationEnvironment, numberOfRepeats, 2, true);
title('Fig. 8 (h)');

