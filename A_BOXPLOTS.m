%boxplots of labels

%load in labels
fileread = 'RGBestsFIX/features_labels.xlsx';

%load in data
R1 = xlsread(fileread,'labels','J2:J18');
ME = xlsread(fileread,'labels','K2:K18');
RE = xlsread(fileread,'labels','L2:L18');
MH = xlsread(fileread,'labels','M2:M18');
RH = xlsread(fileread,'labels','N2:N18');

%boxplot
data = [R1,ME,RE,MH,RH];
figure
boxplot(data)
ylabel('VAS rating')
xtix = {'First rest','Easy maths','Easy rest','Hard maths','Hard rest'}; % Your labels
xtixloc = [1, 2, 3, 4, 5]; % Your label locations
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
print('labels','-depsc')

%two-sample t-tests
[~,pR1ME] = ttest2(R1,ME)
[~,pR1RE] = ttest2(R1,RE)
[~,pR1MH] = ttest2(R1,MH)
[~,pR1RH] = ttest2(R1,RH)
[~,pMERE] = ttest2(ME,RE)
[~,pMEMH] = ttest2(ME,MH)
[~,pMERH] = ttest2(ME,RH)
[~,pREMH] = ttest2(RE,MH)
[~,pRERH] = ttest2(RE,RH)
[~,pMHRH] = ttest2(MH,RH)

