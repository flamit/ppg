%Extracts IBI from raw colour signals imported via Excel

%tidy up
clear
%close all

tic

%load in filename and Ref IBI
fileread = 'RGBestsFIX/5_R5Time2630P7.aviROI_FS_coloursFIX.xlsx';
subject = extractBetween(fileread,'P','.aviROI_FS_coloursFIX.xlsx');
subject = subject{1};
col = char(64+str2double(subject));
IBIfile = 'IBIRefData.xlsx';
IBI_Ref = xlsread(IBIfile,1,strcat(col,':',col));
tfile = 'IBIRefTimes.xlsx';
t_Ref = xlsread(tfile,1,strcat(col,':',col));
%load in sampling frequency
fs = xlsread(fileread,1,'A1');
%sometimes fs is in cell B1
if fs < 10
    fs = xlsread(fileread,1,'B1');
end

%BVP frequency range of interest
LF = 0.7; HF = 4;

%load in raw signals
do_plot = 0;
[red,green,blue] = load_raw(fileread,do_plot);

%calculate the IBIs from the raw colour signals
[t_filtered,IBI_filtered] = raw2ibi(red,green,blue,fs,LF,HF);

%extract features
do_plot = 1;
[AVNN,meanHR,SDNN,RMSSD,pNN50,HRV_HF,HRV_LF,HRV_HFn,HRV_HLn,...
    HRV_LFHFratio] = extract_features(t_filtered,IBI_filtered,do_plot)

export = 0;
if export == 1
    % Get HR data (sliding window)
    [window,HRdata] = slidingHR(t_filtered,IBI_filtered);
    %export to excel file
    filename = strcat('HRdata.xlsx');
    start_position = num2str((str2double(subject)-1)*window + 3);
    xlswrite(filename,HRdata,1,strcat('B',start_position));
    
    % get other data (sliding window)
    do_plot = 0;
    step = 10;
    window = 18;
    features = sliding_features(t_filtered,IBI_filtered,step,window,do_plot);
    %export to excel file
    filename = strcat('Scatter_data.xlsx');
    start_position = num2str((str2double(subject)-1)*window + 3);
    xlswrite(filename,features,1,strcat('B',start_position));
end

%calculate best RMSE of alignment between ref and rPPG
limit = 10;
RMSEs = zeros(1,2*limit+1);

%start with current alignment
RMSE = get_RMSE(IBI_Ref,t_Ref,IBI_filtered,t_filtered);
RMSEs(11) = RMSE;
    
%slide rPPG est back up to ten beats
for step = 1:limit
%trim both sequences
new_ref_IBI = IBI_Ref(1:end-step);
new_ref_t = t_Ref(1:end-step);
new_rPPG_IBI = IBI_filtered(1+step:end);
new_rPPG_t = t_filtered(1+step:end);
RMSE = get_RMSE(new_ref_IBI,new_ref_t,new_rPPG_IBI,new_rPPG_t);
RMSEs(11-step) = RMSE;
end

%slide rPPG est forward up to ten beats
for step = 1:limit
%trim both sequences
new_ref_IBI = IBI_Ref(1+step:end);
new_ref_t = t_Ref(1+step:end);
new_rPPG_IBI = IBI_filtered(1:end-step);
new_rPPG_t = t_filtered(1:end-step);
RMSE = get_RMSE(new_ref_IBI,new_ref_t,new_rPPG_IBI,new_rPPG_t);
RMSEs(11+step) = RMSE;
end

chosen_RMSE = min(RMSEs)

toc



