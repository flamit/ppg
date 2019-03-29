function [window,HRdata] = slidingHR(t_filtered,IBI_filtered)
% Get HR data

step = 15;
window = 18;
HRdata = zeros(window,1);

for w = 1:window
    %selected required interbeats in window
    T = step*(w+1);
    indices1 = logical(t_filtered>T-step);
    t_filtered_2 = t_filtered(indices1);
    indices2 = logical(t_filtered_2<=T+step);
    IBI_window = IBI_filtered(indices1);
    IBI_window = IBI_window(indices2);
    
    %calculate the mean HR and store
    meanIBI = mean(IBI_window);
    HR = 60/meanIBI;
    HRdata(w) = HR;
end