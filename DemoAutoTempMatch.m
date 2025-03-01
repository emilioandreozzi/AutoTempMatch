% DEMO AUTO TEMPLATE MATCHING


% SCG  signal was band-pass filtered in the [7-30] Hz frequency band.
% ECG signal was first band-pass filtered in the [0.5 - 40] Hz frequency
% band. Then, notch filters were used to remove powerline interference and
% its higher harmonics. 

load('data.mat');


%% ECG-FREE HEARTBEATS DETECTION VIA TEMPLATE MATCHING
[Heartbeats,tmpl,wtmpl] = AutoTempMatch(scg,fs);

% VISUALIZATION OF THE SELECTED HEARTBEAT TEMPLATE
normabs = @(x) x/max(abs(x));
figure, plot(t(wtmpl),tmpl,'LineWidth',2), grid on
xlabel('Time [s]'), ylabel('Amplitude [a.u.]')
title('SELECTED TEMPLATE')
axis tight;

% COMPARISON OF HEARTBEATS LOCALIZATION ON ECG AND SCG SIGNALS
figure, plot(t, normabs(ecg)), grid on
hold on
plot(t, normabs(scg))
stem(t(Heartbeats),ones(size(Heartbeats)),'LineStyle','-','Color','k','MarkerFaceColor','r','MarkerEdgeColor','k')
plot(t(Rpeaks), normabs(ecg(Rpeaks)),'.k','MarkerSize',20), 
legend('ECG','SCG','Heartbeats','R-peaks')
xlabel('Time [s]'), ylabel('Amplitude [a.u.]')
title('HEARTBEATS LOCALIZATION')
axis tight;
ylim([-1.2 1.2])



%% COMPARISON WITH ECG ON INTER-BEAT INTERVALS ESTIMATION

% INTER-BEAT INTERVALS
IBIecg = 1000*diff(Rpeaks)/fs;
IBIscg = 1000*diff(Heartbeats)/fs;

% COMPARISON OF INTER-BEAT INTERVALS TRENDS OBTAINED FROM ECG AND SCG
figure
subplot(2,1,1), stem(t(Rpeaks(2:end)), IBIecg,'.b','MarkerSize',20), grid on
xlabel('Time (s)'), ylabel('Inter-beat interval (ms)'), title('Inter-beat intervals from ECG')
subplot(2,1,2), stem(t(Heartbeats(2:end)), IBIscg,'r.','MarkerSize',20), grid on
xlabel('Time (s)'), ylabel('Inter-beat interval (ms)'), title('Inter-beat intervals from SCG')

% COMPARISON OF INTER-BEAT INTERVAL MEASUREMENTS VIA LINEAR REGRESSION,
% CORRELATION, AND BLAND-ALTMAN ANALYSES

% Linear regression and correlation analysis
p = [IBIecg ones(size(IBIecg))]\IBIscg;
slope = p(1);
intercept = p(2);
r = corr(IBIecg,IBIscg);
RLstr = sprintf('y = %1.2f·x + %1.2f\nr = %1.3f',slope,intercept,r);

% Bland-Altman analysis
IBIerr = IBIscg-IBIecg;         % absolute error on IBI estimation
ebias = mean(IBIerr);           % bias: mean error
esd = std(IBIerr);              % standard deviation of errors
LoA = ebias + 1.96*[-esd esd];  % Limits of Agreement (95% confidence interval of errors)
BAstr = sprintf('bias = %1.2f ms\nLoA = [%1.2f %1.2f] ms',ebias,LoA);


% Plots
xa = (min(IBIecg)-20):(max(IBIecg)+20);

figure

% Linear regression plot
subplot(1,2,1)
title('LINEAR REGRESSION PLOT','FontSize',14)
plot(IBIecg,IBIscg,'.','MarkerSize',20), grid on, hold on
plot(xa,xa,'k','LineWidth',1)
plot(xa,slope*xa+intercept,'r--','LineWidth',2)
axis([xa(1) xa(end) xa(1) xa(end)])
xlabel('IBI_{ECG}(ms)','FontSize',20),ylabel('IBI_{SCG}(ms)','FontSize',20)
legend({'Experimental data','Line of equality','Regression line'},'Location','southeast','FontSize',12)
set(gca,'FontSize',14)
text(xa(1)+20,xa(1)+0.94*(xa(end)-xa(1)),RLstr,'FontSize',14)

% Bland-Altman plot
subplot(1,2,2)
title('BLAND-ALTMAN PLOT','FontSize',14)
plot(IBIecg,IBIerr,'.','MarkerSize',20), grid on, hold on
plot(xa,0*xa+ebias,'k','LineWidth',2)
plot(xa,0*xa+LoA(1),'r--',xa,0*xa+LoA(2),'r--','LineWidth',2)
axis([xa(1) xa(end) -40 40])
xlabel('IBI_{ECG}(ms)','FontSize',20),ylabel('ABSOLUTE IBI ERROR (ms)','FontSize',20)
legend({'IBI errors','bias of errors','Limits of Agreement (±1.96 SD)'},'Location','southeast','FontSize',12)
set(gca,'FontSize',14)
text(xa(1)+20,-40+0.94*80,BAstr,'FontSize',14)
