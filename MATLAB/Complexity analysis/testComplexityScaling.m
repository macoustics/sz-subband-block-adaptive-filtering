clear all; close all; clc

RirLength = 700; % [taps]
LoudspeakerFilterLength = 300; % [taps]
NumberOfLoudspeakers = [9];
SourceLowCutoffFrequency = [0];
SourceHighCutoffFrequenc = [24000];
fs = 48000;
NumberOfMicrophones = 2;

NumberOfGradientSteps = 1;
NumberOfLineSearchSteps = 2;

NumberOfChannels = 16;
DecimationFactor = 8;

%% Complexity vs number of microphones
mics = (1:50)*2;
for idx = 1:length(mics)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, mics(idx), NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end

figure
plot(mics,real(fullband) + imag(fullband));
hold on; grid on
plot(mics,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(mics,real(subbandTotal)+imag(subbandTotal));
xlabel('Number of microphones'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')
clear fullband subbandTotal subbandAdaptiveFiltering mics

%% Complexity vs number of loudspeakers
srcs = 2:50;
for idx = 1:length(srcs)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, srcs(idx), NumberOfMicrophones, NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end
figure
plot(srcs,real(fullband) + imag(fullband));
hold on; grid on
plot(srcs,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(srcs,real(subbandTotal)+imag(subbandTotal));
xlabel('Number of loudspeakers'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')

clear fullband subbandTotal subbandAdaptiveFiltering srcs

%% Complexity vs RIR length
rir = 100:100:10000;

for idx = 1:length(rir)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(rir(idx), LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end
figure
plot(rir,real(fullband) + imag(fullband));
hold on; grid on
plot(rir,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(rir,real(subbandTotal)+imag(subbandTotal));
xlabel('RIR length [samples]'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')

clear fullband subbandTotal subbandAdaptiveFiltering rir

%% Complexity vs FIR length
fir = 100:100:10000;

for idx = 1:length(fir)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(RirLength, fir(idx), NumberOfLoudspeakers, NumberOfMicrophones, NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end
figure
plot(fir,real(fullband) + imag(fullband));
hold on; grid on
plot(fir,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(fir,real(subbandTotal)+imag(subbandTotal));
xlabel('FIR length [samples]'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')

clear fullband subbandTotal subbandAdaptiveFiltering fir

%% Complexity vs line search steps
step = 0:100;
for idx = 1:length(step)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), NumberOfGradientSteps, step(idx), SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end
figure
plot(step,real(fullband) + imag(fullband));
hold on; grid on
plot(step,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(step,real(subbandTotal)+imag(subbandTotal));
xlabel('Number of line-search steps'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')

clear fullband subbandTotal subbandAdaptiveFiltering step

%% Complexity vs gradient steps
step = 0:100;
for idx = 1:length(step)
    [tmp1, tmp2, tmp3, tmp4] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, NumberOfChannels, floor(NumberOfChannels*2/4), round(NumberOfChannels*4/(4/2-1)), step(idx), NumberOfGradientSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    fullband(idx) = tmp1;
    subbandAdaptiveFiltering(idx) = tmp2(1);
    subbandTotal(idx) = sum(tmp2);
end
figure
plot(step,real(fullband) + imag(fullband));
hold on; grid on
plot(step,real(subbandAdaptiveFiltering)+imag(subbandAdaptiveFiltering));
plot(step,real(subbandTotal)+imag(subbandTotal));
xlabel('Number of gradient steps'); ylabel('FLOPs')
legend('Fullband','Subband adaptive filtering','Subband total')

clear fullband subbandTotal subbandAdaptiveFiltering step