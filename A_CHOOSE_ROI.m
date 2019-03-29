%Face detection, non-fixed ROI segmentation and BVP extraction

%tidy up
clear
close all
tic
%% UPLOAD VIDEO AND CHOOSE ROI AND BVP FREQUENCY RANGE
%filename
filename = '0807160515609cardio5mins.avi';
ROItrim = 0.2;
LF = 0.7; HF = 4;

%create video object
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

%set up for segmentation
blank = zeros(height,width);

%% FACE DETECTION
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

%% ROI segmentation

%key parameters
seg_height = 10;
seg_width = 10;
cush = 0.1;
limit = 0.5;

preview = 0; %set = 1 to preview ROI
if preview == 1
    figure
end

graph_show = 0; %set = 1 to view PSD of each region

%use whole face to ascertain peak of PSD in range of interest
signal = zeros(1,numFrames);
for k = 1:numFrames
    %get roi for this box segment and frame
    roi = get_roi(allboxPoints,k,blank,ROItrim,...
        1,1,1,1);
    
    %green channel only
    roi = roi(:,:,2);
    roi = uint8(roi);
    imageOI = vidFrames(:,:,2,k).*roi;
    
    %give a preview of roi
    if preview == 1 
        imshow(imageOI)
    end
    
    %get average value of roi (pooling)
    imageOI = imageOI(:);
    imageOI(imageOI==0) = [];
    valueOI = mean(imageOI);
    signal(k) = valueOI;
end
%detrend, normalise and apply bandpass filter to signal
signal = detrend(signal);
signal = (signal-mean(signal))/std(signal);
figure; subplot(4,1,1); plot(signal);
title('Whole face')
%obtain power spectral density estimate using pwelch function
[pxx,f] = pwelch(signal,[],[],[],fs);
subplot(4,1,2); plot(f,pxx)
%plot bandpassed signal and PSD
signal = bandpass(signal,[LF,HF],fs);
subplot(4,1,3); plot(signal);
[pxx,f] = pwelch(signal,[],[],[],fs);
subplot(4,1,4); plot(f,pxx)

%find peak and cushion indices
[p_cent,i_cent] = max(pxx);
f_cent = f(i_cent);
c_low = max(LF,f_cent-cush);
c_high = f_cent+cush; %CHANGE TO MIN()?
%calculate key indices for integration
Lowers_c = find(f>=c_low);
lower_c = Lowers_c(1);
Uppers_c = find(f<=c_high);
upper_c = Uppers_c(end);
Lowers = find(f>=LF);
lower = Lowers(1);
Uppers = find(f<=HF);
upper = Uppers(end);

%set up matrix to store power for segments in range of interest
segment_weights = zeros(seg_height,seg_width);
for hei = 1:seg_height
    for wid = 1:seg_width              
        %run through all frames to collect pulsatility data
        signal = zeros(1,numFrames);
        for k = 1:numFrames
            %get roi for this box segment and frame            
            roi = get_roi(allboxPoints,k,blank,ROItrim,...
                hei,wid,seg_height,seg_width);            
            
            %green channel only
            roi = roi(:,:,2);
            roi = uint8(roi);
            imageOI = vidFrames(:,:,2,k).*roi;
            
            %get average value of roi (pooling)
            imageOI = imageOI(:);
            imageOI(imageOI==0) = [];
            valueOI = mean(imageOI);
            signal(k) = valueOI;
        end
        %detrend and normalise signal
        signal = detrend(signal);
        signal = (signal-mean(signal))/std(signal);
        %plot
        if graph_show == 1
            figure; subplot(4,1,1); plot(signal);
            title(sprintf('Segment height: %d, width: %d',hei,wid))
        end
        %obtain power spectral density estimate using pwelch function
        [pxx,f] = pwelch(signal,[],[],[],fs);
        %plot
        if graph_show == 1
            subplot(4,1,2); plot(f,pxx)
        end
        %bandpassed signal
        signal = bandpass(signal,[LF,HF],fs);
        %plot
        if graph_show == 1
            subplot(4,1,3); plot(signal);
        end
        %PSD
        [pxx,f] = pwelch(signal,[],[],[],fs);
        %plot
        if graph_show == 1
            subplot(4,1,4); plot(f,pxx)
        end
        
        
        %REPLACE WITH get_goodness?
        %key frequencies and power sequences      
        Lowy = pxx(lower:lower_c);
        Cushy = pxx(lower_c:upper_c);
        Highy = pxx(upper_c:upper);
        f_low = f(lower:lower_c);
        f_cush = f(lower_c:upper_c);
        f_high = f(upper_c:upper);       
        
        %integrate to determine goodness metric weight
        int_low = trapz(f_low,Lowy);
        int_cush = trapz(f_cush,Cushy);
        int_high = trapz(f_high,Highy);
        
        goodness = int_cush/(int_low+int_high);
        
        %apply weighting according to goodness of fit
        segment_weights(hei,wid) = goodness;

    end
end

figure
heatmap(segment_weights);
print('ROIseg','-depsc')

toc

