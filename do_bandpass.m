function firBVP = do_bandpass(fs,ALF,AHF,avgBVP,do_plot)
% apply custom bandpass filter (Hamming window) to smoothed source signal

ham_order = 128; %for hamming window size

%fir bandpass filter (Hamming window)
bpFilt = designfilt('bandpassfir','FilterOrder',ham_order,...
    'Window','hamming','CutoffFrequency1',ALF,'CutoffFrequency2',AHF,...
        'SampleRate',fs);
firBVP = filter(bpFilt,avgBVP);

%plot
if do_plot == 1
    figure
    subplot(2,1,1); plot(firBVP)
    title('BVP estimate (after smoothing and bandpassing)')
end

end