function data = scatter_plot(fileread,subject,colour,window)
%creates a scatterplot of HR data

marker_size = 40;
marker = 'o';

data = xlsread(fileread,2,strcat('B',num2str((subject-1)*window + 3),...
    ':C',num2str((subject)*window + 2)));
scatter(data(:,1),data(:,2),marker_size,colour,'filled',marker)

end