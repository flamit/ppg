%Face detection and BVP extraction - using POWER to select source signal
%smooths selected signal first before applying bandpass filter
%uses Hamming window for bandpass filter

%tidy up
clear
close all

tic

%% UPLOAD VIDEO AND CHOOSE ROI AND BVP FREQUENCY RANGE
%filenames
filename = 'Trial1_S1.avi';

%Region of interest
ROIregion = 3; %1=trimmed, 2=eyeless&mouthless, 3=forehead&cheeks, 4=forehead 
ROItrim = 0.2; %proportion of width trimmed from each side

%BVP frequency range of interest
LF = 0.7; HF = 4; 

%create video object and extract frame rate
obj = VideoReader(filename);
fs = obj.FrameRate;

%get the number of frames
numFrames = get(obj, 'numberOfFrames');
frameno = numFrames - 100;

%read all video frames
vidFrames = read(obj);
[height,width,~,~] = size(vidFrames);

%print out three channels
figure
subplot(1,3,1); imshow(vidFrames(:,:,1,frameno));
subplot(1,3,2); imshow(vidFrames(:,:,2,frameno));
subplot(1,3,3); imshow(vidFrames(:,:,3,frameno));

%% FACE DETECTION - Matlab code

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();

% Read a video frame and run the face detector.
videoFileReader = vision.VideoFileReader(filename);
videoFrame      = step(videoFileReader);
bbox            = step(faceDetector, videoFrame);

% Draw the returned bounding box around the detected face.
videoFrame = insertShape(videoFrame, 'Rectangle', bbox);
figure; imshow(videoFrame); title('Detected face');

% Convert the first box into a list of 4 points
% This is needed to be able to visualize the rotation of the object.
bboxPoints = bbox2points(bbox(1, :));

% Detect feature points in the face region.
points = detectMinEigenFeatures(rgb2gray(videoFrame), 'ROI', bbox);

% Display the detected points.
figure, imshow(videoFrame), hold on, title('Detected features');
plot(points);

% Create a point tracker and enable the bidirectional error constraint to
% make it more robust in the presence of noise and clutter.
pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

% Initialize the tracker with the initial point locations and the initial
% video frame.
points = points.Location;
initialize(pointTracker, points, videoFrame);

videoPlayer  = vision.VideoPlayer('Position',...
    [100 100 [size(videoFrame, 2), size(videoFrame, 1)]+30]);

% Make a copy of the points to be used for computing the geometric
% transformation between the points in the previous and the current frames
allboxPoints = zeros(4,2*numFrames);
allboxPoints(:,1:2) = bboxPoints;
oldPoints = points;
frameOI = 1;

while ~isDone(videoFileReader)
    frameOI = frameOI+1;
    % get the next frame
    videoFrame = step(videoFileReader);

    % Track the points. Note that some points may be lost.
    [points, isFound] = step(pointTracker, videoFrame);
    visiblePoints = points(isFound, :);
    oldInliers = oldPoints(isFound, :);
    
    if size(visiblePoints, 1) >= 2 % need at least 2 points
        
        % Estimate the geometric transformation between the old points
        % and the new points and eliminate outliers
        [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
            oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
        
        % Apply the transformation to the bounding box points
        bboxPoints = transformPointsForward(xform, bboxPoints);
        allboxPoints(:,2*frameOI-1:2*frameOI) = bboxPoints;
                
        % Insert a bounding box around the object being tracked
        bboxPolygon = reshape(bboxPoints', 1, []);
        videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, ...
            'LineWidth', 2);
                
        % Display tracked points
        videoFrame = insertMarker(videoFrame, visiblePoints, '+', ...
            'Color', 'white');       
        
        % Reset the points
        oldPoints = visiblePoints;
        setPoints(pointTracker, oldPoints);        
    end
    
    % Display the annotated video frame using the video player object
    step(videoPlayer, videoFrame);
end

% Clean up
release(videoFileReader);
release(videoPlayer);
release(pointTracker);

%% collect average of ROI for red, green and blue channels (POOLING)
preview = 1; %set = 1 to preview ROI
if preview == 1
    figure
end
blank = zeros(height,width);
colours = cell(1,3);
for channel = 1:3
    colour = zeros(1,numFrames);
    for k = 1:numFrames
        %construct ROI using parameters chosen at start
        ROI = construct_ROI(allboxPoints,k,blank,ROItrim,ROIregion);
        
        %one channel only
        ROI = ROI(:,:,2);
        ROI = uint8(ROI);
        imageOI = vidFrames(:,:,channel,k).*ROI;
        
        %give a preview of ROI (green channel)
        if preview == 1 && channel == 2 && k == 1
            imshow(imageOI) 
        end
        
        %get average value of ROI (pooling)
        imageOI = imageOI(:);
        imageOI(imageOI==0) = [];
        valueOI = mean(imageOI);
        colour(k) = valueOI;
    end
    colours{channel} = colour;
end

%% load and preprocess red, green and blue graphs

%raw signals
figure
red = colours{1};
subplot(3,1,1); plot(red, 'r')
title('Raw signals')
green = colours{2};
subplot(3,1,2); plot(green, 'g')
blue = colours{3};
subplot(3,1,3); plot(blue, 'b')

% detrend and normalise revised signals
figure
%detrend
red = detrend(red);
green = detrend(green);
blue = detrend(blue);
%normalise
red = (red-mean(red))/std(red);
green = (green-mean(green))/std(green);
blue = (blue-mean(blue))/std(blue);
%plot
subplot(3,1,1); plot(red, 'r')
title('Detrended and normalised raw signals')
subplot(3,1,2); plot(green, 'g')
subplot(3,1,3); plot(blue, 'b')

%trim to put into report
trimstart = 30;
trimend = 90;
newx1 = ceil(trimstart*fs); newx2 = ceil(trimend*fs);
if newx2<length(red)
    figure
    new_x = (newx1:newx2);
    subplot(3,1,1); plot(new_x,red(newx1:newx2),'r')
    xlabel('frames')
    xlim([newx1,newx2])
    subplot(3,1,2); plot(new_x,green(newx1:newx2),'g')
    xlabel('frames')
    xlim([newx1,newx2])
    subplot(3,1,3); plot(new_x,blue(newx1:newx2),'b')
    xlabel('frames')
    xlim([newx1,newx2])
    print('raw_signals','-depsc')
end

%% compute source signals (using fastICA) - calculate PSD/amp estimates
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
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%TWO CHANNELS: RG
data = [red;green];
start = 4; %start position in matrices above
colour_title = 'RG';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);
%TWO CHANNELS: RB
data = [red;blue];
start = 6; %start position in matrices above
colour_title = 'RB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%TWO CHANNELS: GB
data = [green;blue];
start = 8; %start position in matrices above
colour_title = 'GB';
f_ica = 1; %apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

%ORIGINAL SIGNALS: RGB
data = [red;green;blue];
start = 10; %start position in matrices above
colour_title = 'ORIGINAL red,green,blue';
f_ica = 0; %don't apply fast ICA
[source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica);

[~,selected] = max(pow_peakiness)

%% print out selected PSD and set up custom bandpass filtering

%PSD for selected source signal
noisyBVP = source_signals(selected,:);
[peakiness,centre,band_low,band_high,f,pxx] = PSD_peak(noisyBVP,fs,LF,HF);

% choose max size of narrow band around peak of PSD
bandwidth = 1.6;
%custom bandpass filtering
ALF = max([band_low,LF,centre-0.5*bandwidth]);
AHF = min([HF,band_high,centre+0.5*bandwidth]);

%plot graph
figure
sub1 = subplot(2,1,1); plot(f,pxx)
xlim([0,8])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
sub2 = subplot(2,1,2); plot(f,pxx)
xlim([0.5,4.5])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend_PSD(ALF,AHF,peakiness,centre)

%print separate graphs
hfig = figure;
p_new = copyobj(sub1,hfig);
print('selectedPSD','-depsc')
close
hfig = figure;
p_new = copyobj(sub2,hfig);
legend_PSD(ALF,AHF,peakiness,centre)
print('ZOOMselectedPSD','-depsc')
close

%% smooth selected ('peakiest') source signal 

smoothsteps = 5; %for moving average of source signal

coeffMA = ones(1, smoothsteps)/smoothsteps;
avgBVP = filter(coeffMA, 1, noisyBVP);
figure
subplot(2,1,1); plot(noisyBVP)
hold on
plot(avgBVP)
title('Selected source signal with moving average')

%% apply bandpass filter (Hamming window) to smoothed source signal
ham_order = 128; %for hamming window size

%fir bandpass filter (Hamming window)
bpFilt = designfilt('bandpassfir','FilterOrder',ham_order,...
    'Window','hamming','CutoffFrequency1',ALF,'CutoffFrequency2',AHF,...
        'SampleRate',fs);
firBVP = filter(bpFilt,avgBVP);
subplot(2,1,2); plot(firBVP)
title('BVP estimate (after smoothing and bandpassing)')

%% plot the bandpassed BVP estimate, interpolate, find peaks and calculate IBI

%interpolate signal with a cubic spline function (refines BVP peak point)
sampfr = 256;
samp_freq = 1/sampfr; %sampling frequency 256 Hz
x = (1:length(firBVP))/fs; %time in seconds
xx = 0:samp_freq:length(firBVP)/fs;
interBVP = spline(x,firBVP,xx);
figure
subplot(2,1,1); p1 = plot(x,firBVP);
hold on
p2 = plot(xx,interBVP);
title('Smoothed and interpolated BVP estimate')
legend([p1,p2],{'Smooth','Interpolated'})

%print the interpolated signal for report
opener = ceil(trimstart*sampfr);
closer = ceil(trimend*sampfr);
if closer < length(interBVP)
    sub1 = subplot(2,1,2); plot(xx(opener:closer),interBVP(opener:closer))
    xlim([trimstart,trimend])
    xlabel('time (seconds)')
    %print separate graphs
    hfig = figure;
    p_new = copyobj(sub1,hfig);
    print('BVPest','-depsc')
    close
end

%find peaks of interpolated function with positive and negative orientation
Locations = cell(1,2);
Peaks = zeros(1,2);
figure
for orient = 1:2
    intBVP = ((-1)^orient)*interBVP;
    [~,locs] = findpeaks(intBVP,xx);
    meanCycle = mean(diff(locs));
    peak_tol = 0.65*meanCycle; %min number of frames between each peak
    subplot(2,1,orient); findpeaks(intBVP,xx,'MinPeakDistance',peak_tol);
    if orient == 1
        title('Smoothed and interpolated BVP estimates with peaks')
    end
    [pks,locs] = findpeaks(intBVP,xx,'MinPeakDistance',peak_tol);
    Locations{orient} = locs;
    Peaks(orient) = mean(pks);
end

%choose orientation with highest peaks
[~,indx] = max(Peaks)
locs = Locations{indx};

% calulate and plot IBI
IBI = diff(locs);
t = locs(2:end);
figure
subplot(2,1,1); plot(t,IBI,'.-');
title('Initial IBI of smoothed BVP estimate')

%% apply NC-VT algorithm to filter and remove potential artefacts

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
sub2 = subplot(2,1,2); plot(t_filtered,IBI_filtered,'.-');
xlabel('time (s)')
ylabel('IBIs (s)')
%print separate graphs
hfig = figure;
p_new = copyobj(sub2,hfig);
print('IBI','-depsc')
close
%title('Filtered IBI of smoothed BVP estimate')

%% extract features from IBI

%1. AVNN - average of all NN (IBI) intervals
AVNN = mean(IBI_filtered)

%2. SDNN - standard deviation of all NN (IBI) intervals
SDNN = std(IBI_filtered)

%3. Heart rate estimate
meanHR = 60/AVNN

%obtain power spectral density estimate of IBI using plomb function
[pxx,f] = plomb(IBI_filtered,t_filtered,0.4,'normalized'); %normalised
figure
plot(f,pxx);
xlim([0.04,0.4])
xlabel('Frequency (Hz)')
ylabel('Normalised Power (dB)')
print('IBIPSD','-depsc')
hold on

%calculate key indices for integration
Lowers = find(f>=0.04);
lower = Lowers(1);
Mids = find(f>=0.15);
mid = Mids(1);
Uppers = find(f<=0.4);
upper = Uppers(end);

%key frequencies and power sequences
HFy = pxx(mid:upper);
LFy = pxx(lower:mid);
f_high = f(mid:upper);
f_low = f(lower:mid);

%integrate power to get raw HF and LF values
HRV_HF = trapz(f_high,HFy) %4. unnormalised HF value
HRV_LF = trapz(f_low,LFy) %5. unnormalised LF value
%normalised units
HRV_HFn = HRV_HF*100/(HRV_HF+HRV_LF) %6. normalised HF value
HRV_HLn = HRV_LF*100/(HRV_HF+HRV_LF) %7. normalised LF value

%8. LF/HF ratio
HRV_LFHFratio = HRV_LF/HRV_HF

%peak frequency within HF component
[~,peakHF] = max(HFy);
f_HFpeak = f_high(peakHF);
%plot
location = peakHF + mid - 1;
p3 = plot(f(location),pxx(location),'or');
legend(p3,sprintf('f_{HFpeak} = %.2f',f_HFpeak),'Location','northeast')

%9. estimate of respiration rate from peak frequency within HF component
RR_est = 60*f_HFpeak

%calculate RMSSD and pNN50
IBIdiff = diff(IBI);
RMSSD = sqrt(sum(IBIdiff.^2)/length(IBIdiff)) %10. RMSSD
pNN50 = sum((abs(IBIdiff)>0.05))/length(IBIdiff) %11. pNN50

toc

