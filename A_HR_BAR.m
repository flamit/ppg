%bar chart of HR
close all
clear
fileread = 'GroundTruth_features.xlsx';
P = 16; %number of participants
V = 5; %number of sessions
Ttest = zeros(2*P,V);
All_inputs = zeros(P,V);

for p = 1:P
    start = 10*(p-1) + 3;
    pot_features = xlsread(fileread,'bar',...
        strcat('B',num2str(start),':C',num2str(start+2*V-1)));
    r1 = pot_features(1:2,2);
    R1 = mean(r1);
    %if hard maths task first
    if pot_features(3,1) == 1
        mh = pot_features(3:4,2);
        MH = mean(mh);
        rh = pot_features(5:6,2);
        RH = mean(rh);
        me = pot_features(7:8,2);
        ME = mean(me);
        re = pot_features(9:10,2);
        RE = mean(re);
    %if hard maths task second
    else
        me = pot_features(3:4,2);
        ME = mean(me);
        re = pot_features(5:6,2);
        RE = mean(re);
        mh = pot_features(7:8,2);
        MH = mean(mh);
        rh = pot_features(9:10,2);
        RH = mean(rh);
    end
    
    sessions = [r1,me,re,mh,rh];
    Ttest(2*(p-1)+1:2*p,:) = sessions;
    SESSIONS = [R1,ME,RE,MH,RH];
    All_inputs(p,:) = SESSIONS;
end

%bar graph of actual HR change
do_separated = 0;
if do_separated == 1    
    bar(All_inputs(:,2:5))
    legend('easy maths','easy rest','hard maths','hard rest')
    ylim([0.5,1.5])
    xlabel('Participant no.')
    ylabel('HR/baseline')
    print('HRbar','-depsc')
end

%Mean HRs observed in each session
figure
type = 'RMSSD';
Mean_sesh = mean(Ttest,1);
Std_sesh = std(Ttest,[],1);
bar(Mean_sesh(2:5),'FaceColor','b','EdgeColor','b'); %[1 0.4 0.6]
hold on
%errorbar(Mean_sesh(2:5),Std_sesh(2:5),'r.')
xtix = {'Easy maths','Easy rest','Hard maths','Hard rest'}; % Your labels
xtixloc = [1, 2, 3, 4]; % Your label locations
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc)
xlabel('Session')
ylabel(strcat(type,'/baseline'))
%print('BarmeanRefHRVratio','-depsc')
print(strcat('BarmeanRef',type),'-depsc')

%two-sample t-tests
ME = Ttest(:,2);
RE = Ttest(:,3);
MH = Ttest(:,4);
RH = Ttest(:,5);

[~,pMERE] = ttest2(ME,RE)
[~,pMEMH] = ttest2(ME,MH)
[~,pMERH] = ttest2(ME,RH)
[~,pREMH] = ttest2(RE,MH)
[~,pRERH] = ttest2(RE,RH)
[~,pMHRH] = ttest2(MH,RH)

    
    




