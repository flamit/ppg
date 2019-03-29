function [source_signals,selected] = compute_ICA(red,green,blue,LF,HF,fs,do_plot)
%computes source signals (using fastICA) - calculates PSD estimates
%set up matrices to store source signals and select 'best' one
len = length(red);
source_signals = zeros(12,len);
pow_peakiness = zeros(1,12);

%THREE CHANNELS: RGB
data = [red;green;blue];
start = 1; %start position in matrices above
colour_title = 'RGB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot);

%TWO CHANNELS: RG
data = [red;green];
start = 4; %start position in matrices above
colour_title = 'RG';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot);
%TWO CHANNELS: RB
data = [red;blue];
start = 6; %start position in matrices above
colour_title = 'RB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot);

%TWO CHANNELS: GB
data = [green;blue];
start = 8; %start position in matrices above
colour_title = 'GB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot);

%ORIGINAL SIGNALS: RGB
data = [red;green;blue];
start = 10; %start position in matrices above
colour_title = 'ORIGINAL red,green,blue';
f_ica = 0; %don't apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot);

[~,selected] = max(pow_peakiness)
end