%Imports BVP reference estimate (R5), calculates IBI and extracts features from IBIs
%tidy up
clear
%close all

tic

%% import and plot reference BVP estimates, find peaks and calculate IBI

%key variables
sampfr = 256;
samp_freq = 1/sampfr; %sampling frequency 256 Hz

filereads = {'0726144412P7'};

num_of_vids = length(filereads);
    
for vid_no = 1:num_of_vids
    clearvars -except filereads num_of_vids vid_no sampfr samp_freq
    %load in filename
    fileread = strcat('Ground_Truth/',filereads{vid_no},'.xlsx');
    
    %set up key variables and titles
    R5 = xlsread(fileread,'Inf','P:P')';
    subject = extractBetween(fileread,'P','.xlsx');
    subject = subject{1};
    set_title = strcat('Reference: rest 5, P',subject);
    
    %obtain reference BVP
    refBVP = R5(3:end);
    
    %remove zero entries
    BVPnum = nnz(refBVP);
    refBVP = refBVP(1:BVPnum);
    
    %time axis
    xx = 0:samp_freq:length(refBVP)/sampfr-samp_freq;
    
    %plot
    figure
    subplot(2,1,1); plot(xx,refBVP);
    title(set_title)
    
    %find peaks of reference BVP
    [~,locs] = findpeaks(refBVP,xx);
    %ensure gaps are large enough to avoid mini peaks
    for round = 1:2
        len = ceil(length(locs)/5); %length of top 20% of gaps
        sorted_gaps = sort(diff(locs));
        chosen_gaps = sorted_gaps(end-len:end);
        meanCycle = mean(chosen_gaps);
        peak_tol = 0.65*meanCycle; %min number of frames between each peak
        [~,locs] = findpeaks(refBVP,xx,'MinPeakDistance',peak_tol);
    end
    
    subplot(2,1,2); findpeaks(refBVP,xx,'MinPeakDistance',peak_tol);
    
    % calulate IBI
    IBI = diff(locs);
    t = locs(2:end);
    
    % apply NC-VT algorithm to filter and remove potential artefacts
    thresh_n = 0.2; %step-step threshold
    thresh_m = 0.2; %step-mean threshold
    lenIBI = length(IBI);
    meanIBI = mean(IBI);
    acceptedIBIs = zeros(1,lenIBI);
    %last 'accepted' point
    last = 1;
    accepted = 0;
    for i = 2:lenIBI-1
        if abs(IBI(i)-IBI(last))/IBI(last)<thresh_n ||...
                abs(IBI(i)-IBI(i+1))/IBI(last)<thresh_n||...
                abs(IBI(i)-meanIBI)/meanIBI<thresh_m
            %accept IBI(i)
            last = i;
            accepted = accepted+1;
            acceptedIBIs(accepted) = i;
        end
    end
    
    %remove zero entries
    IBInum = nnz(acceptedIBIs);
    acceptedIBIs = acceptedIBIs(1:IBInum);
    
    %remove rejected IBI values and plot filtered IBI
    IBI_filtered = IBI(acceptedIBIs);
    t_filtered = t(acceptedIBIs);
    figure; subplot(4,1,1);plot(t_filtered,IBI_filtered,'.-'); ylabel('IBIs (s)')
    %title(set_title)
    xlabel('time (s)')
    xlim([50,250])
    ylim([0.7,1.1])
    print('IBIness','-depsc')
    
    %export IBI to excel file
    filename = 'IBIRefData.xlsx';
    xlswrite(filename,IBI_filtered',1,strcat(char(64+str2double(subject)),'1'))
    filename = 'IBIRefTimes.xlsx';
    xlswrite(filename,t_filtered',1,strcat(char(64+str2double(subject)),'1'))
    
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
        xlswrite(filename,{strcat('P',subject)},1,strcat('A',start_position));
        xlswrite(filename,HRdata,1,strcat('C',start_position));
        
        % get other data (sliding window)
        do_plot = 0;
        step = 10;
        window = 18;
        features = sliding_features(t_filtered,IBI_filtered,step,window,do_plot);
        %export to excel file
        filename = strcat('Scatter_data.xlsx');
        start_position = num2str((str2double(subject)-1)*window + 3);
        xlswrite(filename,{strcat('P',subject)},1,strcat('A',start_position));
        xlswrite(filename,features,1,strcat('M',start_position));
    end
    
    toc
end





