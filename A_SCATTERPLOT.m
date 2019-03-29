%Create scatterplots
close all
clear

fileread = 'Scatter_data.xlsx';
window = 18;

%% Side lighting
figure
subjects = {1,2,3,4,22};
colours = {'red','red','blue','red','red'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel('HR using webcam (bpm)')
        ylabel('HR using reference (bpm)')
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print('side','-depsc')

%% Front lighting
figure
subjects = {5,6,7,8,20};
colours = {'red','blue','red','red','red'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel('HR using webcam (bpm)')
        ylabel('HR using reference (bpm)')
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print('front','-depsc')

%% Lamp
figure
subjects = {9,10,18,19,21};
colours = {'blue','blue','red','blue','blue'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel('HR using webcam (bpm)')
        ylabel('HR using reference (bpm)')
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print('lamp','-depsc')

%% Natural light
figure
subjects = {14,15,16,17,23};
colours = {'blue','blue','red','blue','blue'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel('HR using webcam (bpm)')
        ylabel('HR using reference (bpm)')
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print('natural','-depsc')

%% ALL
figure
type = 'HRV LF';
subjects = {1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,20,21,22,23};
colours = {'red','red','red','red','blue','blue','blue','blue','yellow','yellow',...
    'green','green','green','green','yellow','yellow','blue','yellow','red','green'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel(strcat(type,' using webcam (bpm)'))
        ylabel(strcat(type,' using reference (bpm)'))
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print(strcat('all_',type),'-depsc')

%% No lamp
figure
subjects = {1,2,3,4,5,6,7,8,14,15,16,17,20,22,23};
colours = {'red','red','red','red','blue','blue','blue','blue',...
    'green','green','green','green','blue','red','green'};
num_subjects = length(subjects);
DATA = zeros(window*num_subjects,2);

for subj = 1:num_subjects
    subject = subjects{subj};
    colour = colours{subj};
    data = scatter_plot(fileread,subject,colour,window);
    DATA((subj-1)*window+1:window*subj,:) = data;
    if subj == 1
        xlabel('HR using webcam (bpm)')
        ylabel('HR using reference (bpm)')
        hold on
    end
end

%add linear fit and Pearson's cc
lin_fit(DATA)
print('frontnat','-depsc')






