%Test R3 v R5
%tidy up
clear
%close all

tic

%% import and plot reference BVP estimates, find peaks and calculate IBI

%key variables
sampfr = 256;
samp_freq = 1/sampfr; %sampling frequency 256 Hz
files = {'0720202421P1','0725095437P2','0725114340P3','0725135216P4','0726094551P5',...
    '0726114041P6','0726144412P7','0726174523P8','0727094525P9','0727120212P10',...
    '0729122304P14','0729133321P15','0729165929P16','0729185122P17','0730114205P18',...
    '0730133959P19','0802094643P20','0802111708P21','0802131257P22','0802184155P23'};

for subj = 1:length(files)
    file = files{subj};    
    %import reference BVP estimate
    fileread = strcat('Ground_Truth/',file,'.xlsx');
    %get subject number
    subject = extractBetween(fileread,'P','.xlsx');
    subject = subject{1};
    %import data
    R3 = xlsread(fileread,'Inf','J:J')';
    R5 = xlsread(fileread,'Inf','P:P')';
    
    sessions = {R3,R5};
    titles = {'Rest3','Rest5'};
    
    for sesh = 1:length(sessions)
        
        %set graph titles
        set_title = strcat('Reference: ',titles{sesh},', P',subject);
        
        %obtain reference BVP
        refBVP = sessions{sesh}(3:end);
        
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
        
        %plot
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
        figure; plot(t_filtered,IBI_filtered,'.-'); ylabel('IBIs (s)')
        title(set_title)
        
        %extract some features
        AVNN = mean(IBI_filtered);
        meanHR = 60/AVNN
        
        %obtain power spectral density estimate of IBI using plomb function
        [pxx,f] = plomb(IBI_filtered,t_filtered,0.4,'normalized'); %normalised
        figure
        plot(f,pxx);
        xlim([0.04,0.4])
        xlabel('Frequency (Hz)')
        ylabel('Normalised Power (dB)')
        title(set_title)
    end
    
    toc
end


