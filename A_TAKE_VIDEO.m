%record video without autoexposure etc
%installed image acquisition toolbox support package for os generic video interface

%% housekeeping and video set-up
clear
close all

%set up file name
is_test = 0; %set to 1 for tests to avoid clutter
if is_test == 1
    filename = 'test.avi';
else
    %set up timestamp
    format shortg
    c = clock;
    c = fix(c);
    %name of saved file (timestamp will be applied)
    pre_filename = '609cardio5mins';
    filename = sprintf('%02d%02d%02d%02d%02d%s.avi',...
        c(2),c(3),c(4),c(5),c(6),pre_filename);
end

%set up recording time and exposure sensitivity
frames2record = 300; %no of frames to set up frame rate (450f = 15s)
capture = 5; %number of seconds to capture in actual recording (2,400s = 40m)
sensitivity = -6; %sensitivity to light (-10 to -2)

%% run through recording

%create the video input of camera
vid = videoinput('winvideo',1,'MJPG_640x480'); %accessed the laptop's webcam

%get device properties
src = getselectedsource(vid);

%define and set parameters to be changed
properties = {'Exposure','BacklightCompensation','ExposureMode','WhiteBalanceMode'};
values = {sensitivity,'off','manual','manual'};
set(src, properties, values)

%preview video
preview(vid)

%set number of frames to record in preparation shoot
set(vid,'FramesPerTrigger',frames2record);

%record preparation video
start(vid);
wait(vid,Inf);

% retrieve the frames and timestamps for each frame.
[frames,time] = getdata(vid, get(vid,'FramesAvailable'));

% calculate frame rate by averaging difference between each frame's timestamp
framerate = mean(1./diff(time));

%Log frames to disk rather than memory
set(vid,'LoggingMode','disk');

%log to avi file
avi = VideoWriter(filename,'Uncompressed AVI'); %Uncompressed AVI file with RGB24 video
avi.FrameRate = framerate;
set(vid,'DiskLogger',avi);

%set number of frames to record in actual recording
full_frames2record = floor(capture*framerate);
set(vid,'FramesPerTrigger',full_frames2record);

%record full video
start(vid);
wait(vid,Inf);

%housekeeping
video = get(vid,'DiskLogger');
close(video);
delete(vid);
clear vid;
