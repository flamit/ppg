%Extracts features from raw colour signals imported via Excel

%tidy up
clear
close all

tic
%files to run through
filereads = feat_files;
num_of_vids = length(filereads);

%BVP frequency range of interest
LF = 0.7; HF = 4;

for vid_no = 1:num_of_vids
    clearvars -except filereads num_of_vids vid_no LF HF
    %load in filename
    fileread = strcat('RGBestsFIX/',filereads{vid_no},'.xlsx');
    %extract key variables
    subject = extractBetween(fileread,'P','.aviROI_FS_coloursFIX.xlsx');
    subject = subject{1};
    if length(subject)>2
        subject = subject(1:2);
    end
    section = extractBetween(fileread,'FIX/','_');
    section = str2double(section{1});
    try
        difficulty = extractBetween(fileread,'_','Time');
        difficulty = difficulty{1};
    catch
        try
            difficulty = extractBetween(fileread,'_','cutby');
            difficulty = difficulty{1};
        catch
            difficulty = extractBetween(fileread,'_','CUTVERSION');
            difficulty = difficulty{1};
        end
    end
    %load in sampling frequency
    fs = xlsread(fileread,1,'A1');
    %sometimes fs is in cell B1
    if fs < 10
        fs = xlsread(fileread,1,'B1');
    end
    
    %load in raw signals
    do_plot = 0;
    [red,green,blue] = load_raw(fileread,do_plot);
    
    %calculate the IBIs from the raw colour signals
    [t_filtered,IBI_filtered] = raw2ibi(red,green,blue,fs,LF,HF);
    %add a title to the last graph
    if difficulty == 'MH'
        stress_level = 'high stress';
    else
        stress_level = 'low stress';
    end
    title(strcat('P',subject,': section',num2str(section),' -- ',stress_level))
    
    % Extract features (sliding window)
    do_plot = 0;
    step = 60;
    window = 2;
    features = sliding_features(t_filtered,IBI_filtered,step,window,do_plot);
    
    %export data and labels to excel file
    num_sections = 5;
    filename = strcat('RGBestsFIX/features_labels.xlsx');
    %data
    start_position = num2str((str2double(subject)-1)*window*num_sections + (section-1)*window + 3);
    xlswrite(filename,features,'features',strcat('B',start_position));
    %label
    if difficulty == 'MH'
        label = 1; %high stress
    else
        label = 0; % low stress
    end
    xlswrite(filename,[label;label],'features',strcat('M',start_position));
    %participant
    if section == 1
        xlswrite(filename,{strcat('P',subject)},'features',strcat('A',start_position));
    end
end

toc

