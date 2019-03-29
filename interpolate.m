function [xx,interBVP] = interpolate(firBVP,fs,do_plot)
%interpolates signal with a cubic spline function (refines BVP peak point)

sampfr = 256; %sampling frequency 256 Hz
samp_freq = 1/sampfr;

x = (1:length(firBVP))/fs; %time in seconds
xx = 0:samp_freq:length(firBVP)/fs;
interBVP = spline(x,firBVP,xx);

%plot
if do_plot == 1
    figure
    subplot(2,1,1)
    plot(x,firBVP);
    hold on
    plot(xx,interBVP);
    title('Smoothed and interpolated BVP estimate')
end

%print BVP estimate
figure
subplot(4,1,1)
plot(xx,interBVP);
xlim([100,150])
%xlabel('time (s)')
print('BVP_est','-depsc')

end