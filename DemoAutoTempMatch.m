% DEMO AUTO TEMPLATE MATCHING


% SCG  signal was band-pass filtered in the [7-30] Hz frequency band.
% ECG signal was first band-pass filtered in the [0.5 - 40] Hz frequency band. Then, notch filters 
%     were used to remove powerline interference and its higher harmonics.

load('data.mat');


% ECG-free heartbeats detection
[Heartbeats,tmpl,wtmpl] = AutoTempMatch(scg,fs);

% Selected template
normabs = @(x) x/max(abs(x));
figure, plot(t(wtmpl),tmpl,'LineWidth',2), grid on
xlabel('Time [s]'), ylabel('Amplitude [a.u.]')
title('Selected template')
axis tight;

% Comparison with ECG
figure, plot(t, normabs(ecg)), grid on
hold on
plot(t, normabs(scg))
stem(t(Heartbeats),ones(size(Heartbeats)),'LineStyle','-','Color','k','MarkerFaceColor','r','MarkerEdgeColor','k')
plot(t(Rpeaks), normabs(ecg(Rpeaks)),'.k','MarkerSize',20), 
legend('ECG','SCG','Heartbeats','R-peaks')
xlabel('Time [s]'), ylabel('Amplitude [a.u.]')
title('Heartbeats localization')
axis tight;
ylim([-1.2 1.2])
