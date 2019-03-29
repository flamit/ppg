function avgBVP = smooth(noisyBVP,do_plot)
% smooth selected ('peakiest') source signal

smoothsteps = 5; %for moving average of source signal

coeffMA = ones(1, smoothsteps)/smoothsteps;
avgBVP = filter(coeffMA, 1, noisyBVP);

%plot
if do_plot == 1
    figure
    subplot(2,1,1); plot(noisyBVP)
    hold on
    plot(avgBVP)
    title('Selected source signal with moving average')
end

end
