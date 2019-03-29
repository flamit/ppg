function locs = find_peaks(interBVP,xx,do_plot)
%find peaks of interpolated function and choose orientation with highest peaks

%find peaks of interpolated function with positive and negative orientation
Locations = cell(1,2);
Peaks = zeros(1,2);
if do_plot == 1
    figure
end
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
    if do_plot == 1
        subplot(2,1,orient); findpeaks(intBVP,xx,'MinPeakDistance',peak_tol);
        if orient == 1
            title('Smoothed and interpolated BVP estimates with peaks')
        end
    end
    Locations{orient} = locs;
    Peaks(orient) = mean(pks);
end

%choose orientation with highest peaks
[~,indx] = max(Peaks)
locs = Locations{indx};
