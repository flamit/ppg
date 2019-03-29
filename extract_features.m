function [AVNN,meanHR,SDNN,RMSSD,pNN50,HRV_HF,HRV_LF,HRV_HFn,HRV_HLn,...
    HRV_LFHFratio] = extract_features(t_filtered,IBI_filtered,do_plot)
% extracts features from IBI

%1. AVNN - average of all NN (IBI) intervals
AVNN = mean(IBI_filtered);

%2. Heart rate estimate
meanHR = 60/AVNN;

%3. SDNN - standard deviation of all NN (IBI) intervals
SDNN = std(IBI_filtered);

%calculate RMSSD and pNN50
IBIdiff = diff(IBI_filtered);
RMSSD = sqrt(sum(IBIdiff.^2)/length(IBIdiff)); %4. RMSSD
pNN50 = sum((abs(IBIdiff)>0.05))/length(IBIdiff); %5. pNN50

%obtain power spectral density estimate of IBI using plomb function
[pxx,f] = plomb(IBI_filtered,t_filtered,0.4,'normalized'); %normalised

%plot
if do_plot == 1
    figure
    plot(f,pxx);
    xlim([0.04,0.4])
    xlabel('Frequency (Hz)')
    ylabel('Normalised Power (dB)')
    ylim([0,18])
    print('Plomb','-depsc')
end

%calculate key indices for integration
Lowers = find(f>=0.04);
lower = Lowers(1);
Mids = find(f>=0.15);
mid = Mids(1);
Uppers = find(f<=0.4);
upper = Uppers(end);

%key frequencies and power sequences
HFy = pxx(mid:upper);
LFy = pxx(lower:mid);
f_high = f(mid:upper);
f_low = f(lower:mid);

%integrate power to get raw HF and LF values
HRV_HF = trapz(f_high,HFy); %6. unnormalised HF value
HRV_LF = trapz(f_low,LFy); %6. unnormalised LF value
%normalised units
HRV_HFn = HRV_HF*100/(HRV_HF+HRV_LF); %8. normalised HF value
HRV_HLn = HRV_LF*100/(HRV_HF+HRV_LF); %9. normalised LF value

%10. LF/HF ratio
HRV_LFHFratio = HRV_LF/HRV_HF;

end
