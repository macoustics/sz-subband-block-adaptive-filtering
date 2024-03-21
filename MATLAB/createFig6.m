clear all; close all; clc

RirLength = 700; % [taps]
LoudspeakerFilterLength = 300; % [taps]
NumberOfLoudspeakers = [9;11;7];
SourceLowCutoffFrequency = [100; 1500; 4000];
SourceHighCutoffFrequenc = [1500; 4000; 18000];
fs = 48000;
NumberOfMicrophones = 6; % 3 mics in the bright and dark zones respectively

NumberOfGradientSteps = 1;
NumberOfLineSearchSteps = 2;

ChannelNumberList = 2.^(2:7);

CostFullband = zeros(length(ChannelNumberList),1);
CostSubband = zeros(length(ChannelNumberList),1);
CostSubbandActiveChannel = zeros(length(ChannelNumberList),1);

for channelIdx = 1:length(ChannelNumberList)
    [tmp1, tmp2, tmp3] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, ChannelNumberList(channelIdx), round(ChannelNumberList(channelIdx)/2), ChannelNumberList(channelIdx)*4, NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    CostFullband(channelIdx) = sum(tmp1); 
    CostSubband1(channelIdx) = sum(tmp2);
    CostSubbandActiveChannel1(channelIdx) = sum(tmp3);
    CostSBAC1_AF(channelIdx) = tmp3(1);
    CostSBAC1_An(channelIdx) = tmp3(2);
    CostSBAC1_Sn(channelIdx) = tmp3(3);
    [tmp1, tmp2, tmp3] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, ChannelNumberList(channelIdx), round(ChannelNumberList(channelIdx)/1.5), ChannelNumberList(channelIdx)*8, NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    CostSubband2(channelIdx) = sum(tmp2);
    CostSubbandActiveChannel2(channelIdx) = sum(tmp3);
    CostSBAC2_AF(channelIdx) = tmp3(1);
    CostSBAC2_An(channelIdx) = tmp3(2);
    CostSBAC2_Sn(channelIdx) = tmp3(3);
    [tmp1, tmp2, tmp3] = EvaluateSubbandComplexity(RirLength, LoudspeakerFilterLength, NumberOfLoudspeakers, NumberOfMicrophones, ChannelNumberList(channelIdx), round(ChannelNumberList(channelIdx)*1/(1.25)), ChannelNumberList(channelIdx)*16, NumberOfGradientSteps, NumberOfLineSearchSteps, SourceLowCutoffFrequency, SourceHighCutoffFrequenc, fs);
    CostSubband3(channelIdx) = sum(tmp2);
    CostSubbandActiveChannel3(channelIdx) = sum(tmp3);
    CostSBAC3_AF(channelIdx) = tmp3(1);
    CostSBAC3_An(channelIdx) = tmp3(2);
    CostSBAC3_Sn(channelIdx) = tmp3(3);
end

figure
plot(2:7,real(CostFullband)+imag(CostFullband),'k','linewidth',1);
hold on; grid on;
plot(2:7,real(CostSubband1)+imag(CostSubband1),'b','linewidth',1)
plot(2:7,real(CostSubbandActiveChannel1)+imag(CostSubbandActiveChannel1),'r','linewidth',1)
plot(2:7,real(CostSubband2)+imag(CostSubband2),'--b','linewidth',1)
plot(2:7,real(CostSubbandActiveChannel2)+imag(CostSubbandActiveChannel2),'--r','linewidth',1)
plot(2:7,real(CostSubband3)+imag(CostSubband3),'-.b','linewidth',1)
plot(2:7,real(CostSubbandActiveChannel3)+imag(CostSubbandActiveChannel3),'-.r','linewidth',1)
legend('No subbands','All bands: L_p=4K, D=K/2','Only active bands: L_p=4K, R=2','All bands: L_p=8K, R=3/2', 'Only active bands bands: L_p=8K, R=3/2', 'All bands: L_p=16K, R=5/4', 'Only active bands: L_p=16K, R=5/4')
xlabel('Number of filter bank channels [2^x]')
ylabel('Complexity in flops per output sample')
set(gca,'yscale','log')

