
%% Contrast Plot
openfig('ContrastSSR.fig')
title('')
legend off
set(gca,'FontSize',8)
% ylabel('Contrast [dB]'); xlabel('Frequency [Hz]')
% ylim([-20 2]); xlim([0 1172])

si = [7.5,5];
% si = [7.5,4.5];

set(gcf,'units','centimeters');
set(gca,'units','centimeters');
pF = get(gcf,'Position');

% Set plot dimensions (and add some margin for labels + cbar)
set(gca,'Position',[0, 0, si(1)-1.5, si(2)-1.5])

% Evaluate the TightInset margins (changes every time the plot is scaled)
ti = get(gca, 'TightInset');
% Get the current dimensions of the plot
pA = get(gca,'Position');
% Adjust the plot to the left / bottom TightInset margins
set(gca,'Position',[ti(1), ti(2), pA(3), pA(4)])
% Adjust the dimensions of the figure size
set(gcf,'Position',[pF(1),pF(2),si(1),si(2)]);
set(gca, 'XTick', [20 50 100 200]);
set(gca, 'YTick', [-10 0 10 20 30 40])

% print('ContrastSSR','-depsc');

%% SVD plot
openfig('SingularValuesSSR.fig')
title('')
set(gca,'FontSize',8)
ylabel('\sigma_l / \sigma_1'); xlabel('Frequency [Hz]')
ylim([1e-3 1e0]); xlim([20 300])

si = [7.5,5.5];

set(gcf,'units','centimeters');
set(gca,'units','centimeters');
pF = get(gcf,'Position');

% Set plot dimensions (and add some margin for labels + cbar)
set(gca,'Position',[0, 0, si(1)-1.5, si(2)-1.5])

% Evaluate the TightInset margins (changes every time the plot is scaled)
ti = get(gca, 'TightInset');
% Get the current dimensions of the plot
pA = get(gca,'Position');
% Adjust the plot to the left / bottom TightInset margins
set(gca,'Position',[ti(1), ti(2), pA(3), pA(4)])
% Adjust the dimensions of the figure size
set(gcf,'Position',[pF(1),pF(2),si(1),si(2)]);
set(gca, 'XTick', [20 50 100 200]);
set(gca, 'YTick', [1e-3 1e-2 1e-1 1e0])

% print('svdSSR','-depsc');