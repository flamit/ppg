function features = sliding_features(t_filtered,IBI_filtered,step,window,do_plot)
% Extract features via a sliding window

features = zeros(window,10);

for w = 1:window
    %selected required interbeats in window
    T = step*(w+1);
    indices1 = logical(t_filtered>T-step);
    t_window = t_filtered(indices1);
    indices2 = logical(t_window<=T+step);
    t_window = t_window(indices2);
    IBI_window = IBI_filtered(indices1);
    IBI_window = IBI_window(indices2);
    
    %extract the features and store
    [AVNN,meanHR,SDNN,RMSSD,pNN50,HRV_HF,HRV_LF,HRV_HFn,HRV_HLn,...
        HRV_LFHFratio] = extract_features(t_window,IBI_window,do_plot);
    
    features(w,1) = AVNN;
    features(w,2) = meanHR;
    features(w,3) = SDNN;
    features(w,4) = RMSSD;
    features(w,5) = pNN50;
    features(w,6) = HRV_HF;
    features(w,7) = HRV_LF;
    features(w,8) = HRV_HFn;
    features(w,9) = HRV_HLn;
    features(w,10) = HRV_LFHFratio;

end

end