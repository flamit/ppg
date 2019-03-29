function [red,green,blue] = load_raw(fileread,do_plot)
%loads in raw signals from an Excel file

%load raw signals
red = xlsread(fileread,1,'C:C')';
green = xlsread(fileread,1,'D:D')';
blue = xlsread(fileread,1,'E:E')';

%plot
if do_plot == 1
    figure
    subplot(3,1,1); plot(red, 'r')
    title('Raw signals')
    subplot(3,1,2); plot(green, 'g')
    subplot(3,1,3); plot(blue, 'b')
end

end