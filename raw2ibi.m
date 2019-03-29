function [t_filtered,IBI_filtered] = raw2ibi(red,green,blue,fs,LF,HF)
%calculate the IBIs from the raw colour signals

% detrend and normalise raw signals
do_plot = 0;
[red,green,blue] = detrend_norm(red,green,blue,do_plot);

% compute source signals (using fastICA) - calculate PSD estimates
do_plot = 0;
[source_signals,selected] = compute_ICA(red,green,blue,LF,HF,fs,do_plot);

% print out selected PSD and set up custom bandpass filtering
do_plot = 1;
[noisyBVP,ALF,AHF] = selected_PSD(LF,HF,fs,source_signals,selected,do_plot);

% smooth selected ('peakiest') source signal 
do_plot = 0;
avgBVP = smooth(noisyBVP,do_plot);

% apply bandpass filter (Hamming window) to smoothed source signal
do_plot = 0;
firBVP = do_bandpass(fs,ALF,AHF,avgBVP,do_plot);

%interpolate signal with a cubic spline function (refines BVP peak point)
do_plot = 1;
[xx,interBVP] = interpolate(firBVP,fs,do_plot);

%find peaks of interpolated function and choose orientation with highest peaks
do_plot = 0;
locs = find_peaks(interBVP,xx,do_plot);

% calulate IBI and apply NC-VT algorithm
do_plot = 1; do_separate = 1;
[t_filtered,IBI_filtered] = IBI_NCVT(locs,do_plot,do_separate);

end