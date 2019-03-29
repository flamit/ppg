function RMSE = get_RMSE(IBI_Ref,t_Ref,IBI_filtered,t_filtered)
%calculates RMSE of alignment between ref and rPPG IBI estimates
reflen = length (IBI_Ref);
squared_diffs = zeros(1,reflen);
for beat = 1:reflen
    %consider time of reference beat
    t = t_Ref(beat);
    refval = IBI_Ref(beat);
    %look at closet beat to this in rPPG example
    [~,index] = min(abs(t_filtered-t));
    rPPGval = IBI_filtered(index);
    %calculate squared difference
    squared_diff = (refval-rPPGval)^2;
    squared_diffs(beat) = squared_diff;
end
%calculate RMSE
RMSE = sqrt(sum(squared_diffs)/reflen);
end