function [red,green,blue] = detrend_norm(red,green,blue,do_plot)
% detrends and normalises raw signals

%detrend
red = detrend(red);
green = detrend(green);
blue = detrend(blue);
%normalise
red = (red-mean(red))/std(red);
green = (green-mean(green))/std(green);
blue = (blue-mean(blue))/std(blue);
%plot
if do_plot == 1
    figure
    subplot(3,1,1); plot(red, 'r')
    title('Detrended and normalised raw signals')
    subplot(3,1,2); plot(green, 'g')
    subplot(3,1,3); plot(blue, 'b')
end
end