function add_legend(ALF,AHF,peakiness,centre)
%Adds a legend to the zoomed PSD plot
hold on; p1 = plot([ALF,ALF],[0,peakiness],'--r'); plot([AHF,AHF],[0,peakiness],'--r');
p2 = plot(centre,0,'xg'); hold off
legend([p2,p1],'peak HR estimate','custom bandwidth')
end
