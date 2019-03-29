function [peakiness,centre,band_low,band_high,f,pxx,adLB,adUB] = PSD_peak(source,fs,LF,HF)
% Measure peakiness of power in freq. range of interest

%obtain power spectral density estimate using pwelch function
[pxx,f] = pwelch(source,[],[],[],fs);

%get lower bound and upper bound indices for frequency range of interest
LBs = find(f<=LF);
LB = LBs(end);
UBs = find(f>=HF);
UB = UBs(1);

%measure 'peakiness' of power within frequency range of interest
findpeaks(pxx,f)
xlim([0.7,4.0])
poi = pxx(LB:UB);
[pks,locs] = findpeaks(poi);
[peakiness,place] = max(pks);
%find centre of peak
point = locs(place);
centre_indx = LB+point-1;
centre = f(centre_indx);

%ADAPTIVE BANDPASS FILTER SETUP
span = 2;
no_peaks = 6;

while span > 1.6 && no_peaks > 1
    %consider top five or so peaks
    no_peaks = no_peaks - 1;
    %get the indices which are highest
    [~, sortIndex] = sort(pks(:), 'descend');  %sort the values in descending order
    maxIndices = sortIndex(1:no_peaks);  %index of the largest values
    maxLocations = LB - 1 + locs(maxIndices);
    maxFrequencies = f(maxLocations)
    low_freq = min(maxFrequencies);
    high_freq = max(maxFrequencies);
    adLB = max(LF,low_freq-0.3);
    adUB = min(HF,high_freq+0.3);
    span = adUB - adLB
end
    
%CUSTOM BANDPASS FILTER SETUP
%calculate adaptive bandwidth to consider around peak
pk_threshold = 0.4*peakiness;
drop_threshold = 0.2*peakiness;

%lowest peak to consider
low_pk_i = place;
while low_pk_i > 1 && pks(low_pk_i - 1) > pk_threshold 
    low_pk_i = low_pk_i - 1;
end
point = locs(low_pk_i); 
low_pk_indx = LB+point-1;
%at which frequency does that drop down to a low value PSD
poss_lowers = find(pxx(1:low_pk_indx)<drop_threshold); 
lower = max(poss_lowers);
band_low = f(lower);

%highest peak to consider
high_pk_i = place;
while high_pk_i < length(pks) && pks(high_pk_i + 1) > pk_threshold 
    high_pk_i = high_pk_i + 1;
end
point = locs(high_pk_i); 
high_pk_indx = LB+point-1;
%at which frequency does that drop down to a low value PSD
poss_highers = find(pxx(high_pk_indx:end)<drop_threshold); 
higher = high_pk_indx + min(poss_highers) - 1;
band_high = f(higher);

end