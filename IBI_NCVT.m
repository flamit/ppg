function [t_filtered,IBI_filtered] = IBI_NCVT(locs,do_plot,do_separate)
% calulates IBI and applies NC-VT algorithm

%calculate intial IBI
IBI = diff(locs);
t = locs(2:end);
if do_plot == 1
    figure
    subplot(2,1,1); plot(t,IBI,'.-');
    title('Initial IBI of smoothed BVP estimate')
end

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
if do_plot == 1
    subplot(2,1,2); plot(t_filtered,IBI_filtered,'.-');
    xlabel('time (s)')
    ylabel('IBIs (s)')
end

%plot separate graphs
if do_separate == 1
    figure; subplot(4,1,1);
    plot(t_filtered,IBI_filtered,'.-');
    %xlabel('time (s)')
    xlim([50,250])
    ylim([0.7,1.1])
    ylabel('IBIs (s)')
    print('IBIness','-depsc')
    %title('Filtered IBI of smoothed BVP estimate')
end

end