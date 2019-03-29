%Extracts IBI from raw colour signals imported via Excel

%tidy up
clear
%close all

tic

%% load and preprocess red, green and blue graphs

%load in filename
filename = 'RGBestsFIX/5_R5Time2630P7.aviROI_FS_coloursFIX.xlsx';

%BVP frequency range of interest
LF = 0.7; HF = 4;

%load in sampling frequency
fs = xlsread(filename,1,'A1');

%sometimes fs is in cell B1
if fs < 10
    fs = xlsread(filename,1,'B1');
end

%load in raw signals
figure
red = xlsread(filename,1,'C:C')';
subplot(3,1,1); plot(red, 'r')
title('Raw signals')
green = xlsread(filename,1,'D:D')';
subplot(3,1,2); plot(green, 'g')
blue = xlsread(filename,1,'E:E')';
subplot(3,1,3); plot(blue, 'b')

% detrend and normalise revised signals
figure
%detrend
red = detrend(red);
green = detrend(green);
blue = detrend(blue);
%normalise
red = (red-mean(red))/std(red);
green = (green-mean(green))/std(green);
blue = (blue-mean(blue))/std(blue);
%plot
subplot(3,1,1); plot(red, 'r')
title('Detrended and normalised raw signals')
subplot(3,1,2); plot(green, 'g')
subplot(3,1,3); plot(blue, 'b')

%% compute source signals (using fastICA) - calculate PSD/amp estimates
%set up matrices to store source signals and select 'best' one
len = length(red);
source_signals = zeros(12,len);
pow_peakiness = zeros(1,12);

%THREE CHANNELS: RGB
data = [red;green;blue];
start = 1; %start position in matrices above
colour_title = 'RGB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%TWO CHANNELS: RG
data = [red;green];
start = 4; %start position in matrices above
colour_title = 'RG';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);
%TWO CHANNELS: RB
data = [red;blue];
start = 6; %start position in matrices above
colour_title = 'RB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%TWO CHANNELS: GB
data = [green;blue];
start = 8; %start position in matrices above
colour_title = 'GB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%ORIGINAL SIGNALS: RGB
data = [red;green;blue];
start = 10; %start position in matrices above
colour_title = 'ORIGINAL red,green,blue';
f_ica = 0; %don't apply fast ICA
print_out = 0;
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

[~,selected] = max(pow_peakiness)

%% print out selected PSD and set up custom bandpass filtering

%PSD for selected source signal
noisyBVP = source_signals(selected,:);
[peakiness,centre,band_low,band_high,f,pxx] = PSD_peak(noisyBVP,fs,LF,HF);

% choose max size of narrow band around peak of PSD
bandwidth = 1.6;
%custom bandpass filtering
ALF = max([band_low,LF,centre-0.5*bandwidth]);
AHF = min([HF,band_high,centre+0.5*bandwidth]);

%plot graph
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

%print separate graphs
hfig = figure;
p_new = copyobj(sub1,hfig);
print('selectedPSD','-depsc')
close
hfig = figure;
p_new = copyobj(sub2,hfig);
legend_PSD(ALF,AHF,peakiness,centre)
print('ZOOMselectedPSD','-depsc')
close

%% smooth selected ('peakiest') source signal 

smoothsteps = 5; %for moving average of source signal

coeffMA = ones(1, smoothsteps)/smoothsteps;
avgBVP = filter(coeffMA, 1, noisyBVP);
figure
subplot(2,1,1); plot(noisyBVP)
hold on
plot(avgBVP)
title('Selected source signal with moving average')

%% calculate BVP and IBI from bandpassed signal

% apply bandpass filter (Hamming window) to smoothed source signal
ham_order = 128; %for hamming window size

%fir bandpass filter (Hamming window)
bpFilt = designfilt('bandpassfir','FilterOrder',ham_order,...
    'Window','hamming','CutoffFrequency1',ALF,'CutoffFrequency2',AHF,...
        'SampleRate',fs);
firBVP = filter(bpFilt,avgBVP);
figure
subplot(2,1,1); plot(firBVP)
title('BVP estimate (after smoothing and bandpassing)')

%interpolate signal with a cubic spline function (refines BVP peak point)
sampfr = 256;
samp_freq = 1/sampfr; %sampling frequency 256 Hz
x = (1:length(firBVP))/fs; %time in seconds
xx = 0:samp_freq:length(firBVP)/fs;
interBVP = spline(x,firBVP,xx);
figure
p1 = plot(x,firBVP);
hold on
p2 = plot(xx,interBVP);
title('Smoothed and interpolated BVP estimate')
legend([p1,p2],{'Smooth','Interpolated'})

%find peaks of interpolated function with positive and negative orientation
Locations = cell(1,2);
Peaks = zeros(1,2);
figure
for orient = 1:2
    intBVP = ((-1)^orient)*interBVP;
    [~,locs] = findpeaks(intBVP,xx);
    %ensure gaps are large enough to avoid mini peaks
    for round = 1:2
        len = ceil(length(locs)/5); %length of top 20% of gaps
        sorted_gaps = sort(diff(locs));
        chosen_gaps = sorted_gaps(end-len:end);
        meanCycle = mean(chosen_gaps);
        peak_tol = 0.65*meanCycle; %min number of frames between each peak
        [pks,locs] = findpeaks(intBVP,xx,'MinPeakDistance',peak_tol);
    end
    %plot
    subplot(2,1,orient); findpeaks(intBVP,xx,'MinPeakDistance',peak_tol);
    if orient == 1
        title('Smoothed and interpolated BVP estimates with peaks')
    end
    Locations{orient} = locs;
    Peaks(orient) = mean(pks);
end

%choose orientation with highest peaks
[~,indx] = max(Peaks)
locs = Locations{indx};

% calulate and plot IBI
IBI = diff(locs);
t = locs(2:end);
figure
subplot(2,1,1); plot(t,IBI,'.-');
title('Initial IBI of smoothed BVP estimate')

% apply NC-VT algorithm to filter and remove potential artefacts

thresh_n = 0.2; %step-step threshold
thresh_m = 0.2; %step-mean threshold
lenIBI = length(IBI);
meanIBI = mean(IBI);
acceptedIBIs = zeros(1,lenIBI);
%last 'accepted' point
last = 1;
accepted = 0;
for i = 2:lenIBI-1
    if abs(IBI(i)-IBI(last))/IBI(last)<thresh_n ||...
            abs(IBI(i)-IBI(i+1))/IBI(last)<thresh_n||...
            abs(IBI(i)-meanIBI)/meanIBI<thresh_m
        %accept IBI(i)
        last = i;
        accepted = accepted+1;
        acceptedIBIs(accepted) = i;
    end
end

%remove zero entries
IBInum = nnz(acceptedIBIs);
acceptedIBIs = acceptedIBIs(1:IBInum);

%remove rejected IBI values and plot filtered IBI
IBI_filtered = IBI(acceptedIBIs);
t_filtered = t(acceptedIBIs);
subplot(2,1,2); plot(t_filtered,IBI_filtered,'.-');
xlabel('time (s)')
ylabel('IBIs (s)')

%print separate graphs
figure
plot(t_filtered,IBI_filtered,'.-');
xlabel('time (s)')
ylabel('IBIs (s)')
title('Filtered IBI of smoothed BVP estimate')

%% extract features from IBI

%1. AVNN - average of all NN (IBI) intervals
AVNN = mean(IBI_filtered)

%2. SDNN - standard deviation of all NN (IBI) intervals
SDNN = std(IBI_filtered)

%3. Heart rate estimate
meanHR = 60/AVNN

%obtain power spectral density estimate of IBI using plomb function
[pxx,f] = plomb(IBI_filtered,t_filtered,0.4,'normalized'); %normalised
figure
plot(f,pxx);
xlim([0.04,0.4])
xlabel('Frequency (Hz)')
ylabel('Normalised Power (dB)')
hold on

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
HRV_HF = trapz(f_high,HFy) %4. unnormalised HF value
HRV_LF = trapz(f_low,LFy) %5. unnormalised LF value
%normalised units
HRV_HFn = HRV_HF*100/(HRV_HF+HRV_LF) %6. normalised HF value
HRV_HLn = HRV_LF*100/(HRV_HF+HRV_LF) %7. normalised LF value

%8. LF/HF ratio
HRV_LFHFratio = HRV_LF/HRV_HF

%peak frequency within HF component
[~,peakHF] = max(HFy);
f_HFpeak = f_high(peakHF);
%plot
location = peakHF + mid - 1;
p3 = plot(f(location),pxx(location),'or');
legend(p3,sprintf('f_{HFpeak} = %.2f',f_HFpeak),'Location','northwest')

%9. estimate of respiration rate from peak frequency within HF component
RR_est = 60*f_HFpeak

%calculate RMSSD and pNN50
IBIdiff = diff(IBI);
RMSSD = sqrt(sum(IBIdiff.^2)/length(IBIdiff)) %10. RMSSD
pNN50 = sum((abs(IBIdiff)>0.05))/length(IBIdiff) %11. pNN50

toc

