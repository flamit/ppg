%Imports BVP estimate, calculates IBI and extracts and exports features from IBIs
%tidy up
clear
close all

tic

%% import and plot reference BVP estimates, find peaks and calculate IBI

%key variables
sampfr = 256;
samp_freq = 1/sampfr; %sampling frequency 256 Hz
num_sections = 5;
num_features = 10;
step = 60;
window = 2;
hard_first = [2,4,6,8,10,14,16,18,20,22];

%import reference BVP estimate
fileread = 'Ground_Truth/0727120212P10.xlsx';
subject = extractBetween(fileread,'P','.xlsx');
subject = subject{1};
R1 = xlsread(fileread,'Inf','D:D')';
M2 = xlsread(fileread,'Inf','G:G')';
R3 = xlsread(fileread,'Inf','J:J')';
M4 = xlsread(fileread,'Inf','M:M')';
R5 = xlsread(fileread,'Inf','P:P')';

sessions = {R1,M2,R3,M4,R5};
titles = {'Rest1','Maths2','Rest3','Maths4','Rest5'};
Full_features = zeros(num_sections*2,num_features);

for sesh = 1:num_sections
    
    %obtain reference BVP
    refBVP = sessions{sesh}(3:end);
    
    %remove zero entries
    BVPnum = nnz(refBVP);
    refBVP = refBVP(1:BVPnum);
    
    %time axis
    xx = 0:samp_freq:length(refBVP)/sampfr-samp_freq;
    
    %plot
    figure
    subplot(3,1,1); plot(xx,refBVP);
    title(titles{sesh})
    
    %find peaks of interpolated function and choose orientation with highest peaks
    do_plot = 0;
    locs = find_peaks(refBVP,xx,do_plot);
    
    % calulate IBI and apply NC-VT algorithm
    do_plot = 0; do_separate = 1;
    [t_filtered,IBI_filtered] = IBI_NCVT(locs,do_plot,do_separate);
    
    % Extract features (sliding window)
    do_plot = 0;    
    features = sliding_features(t_filtered,IBI_filtered,step,window,do_plot);
    
    Full_features(2*(sesh-1)+1:2*sesh,:) = features;
end

%export data and labels to excel file
filename = strcat('GroundTruth_features.xlsx');
%data
start_position = num2str((str2double(subject)-1)*window*num_sections + 3);
xlswrite(filename,Full_features,'features',strcat('B',start_position));
%label
if any(hard_first == str2double(subject))
    place = 3; %high stress
else
    place = 7; % low stress
end
label = zeros(window*num_sections,1);
label([place,place+1]) = 1;
xlswrite(filename,label,'features',strcat('M',start_position));
%participant
xlswrite(filename,{strcat('P',subject)},'features',strcat('A',start_position));


toc


