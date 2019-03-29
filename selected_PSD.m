function [noisyBVP,ALF,AHF] = selected_PSD(LF,HF,fs,source_signals,selected,do_plot)
% print out selected PSD and set up custom bandpass filtering

%PSD for selected source signal
noisyBVP = source_signals(selected,:);
[peakiness,centre,band_low,band_high,f,pxx,adLB,adUB] = PSD_peak(noisyBVP,fs,LF,HF);

%CHOOSE FILTER TYPE
filter_type = 1; %1 = custom, 2 = adaptive, 3 = full

if filter_type == 1 %custom
% choose max size of narrow band around peak of PSD
bandwidth = 1.6;
%custom bandpass filtering
ALF = max([band_low,LF,centre-0.5*bandwidth]);
AHF = min([HF,band_high,centre+0.5*bandwidth]);
elseif filter_type == 2 %adaptive
    ALF = adLB;
    AHF = adUB;    
else %full
    ALF = LF;
    AHF = HF;
end

%plot graph
if do_plot == 1
    figure
    sub1 = subplot(2,1,1); plot(f,pxx)
    xlim([0,8])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    sub2 = subplot(2,1,2); plot(f,pxx)
    xlim([0.5,4.5])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    legend_PSD(ALF,AHF,peakiness,centre)
end

%print PSD plot
figure
subplot(4,1,1); plot(f,pxx)
xlim([0.5,4.5])
%xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend_PSD(ALF,AHF,peakiness,centre)
print('PSDness','-depsc')

end

