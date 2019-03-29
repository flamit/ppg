function lin_fit(DATA)
%adds a linear fit and Pearson's cc

%line fit
p = polyfit(DATA(:,1),DATA(:,2),1);
x1 = min(DATA(:,1)); x2 = max(DATA(:,1));
x = linspace(x1,x2);
y = polyval(p,x);
plot(x,y,'k-.','LineWidth',1)
%add text to plot indicating Pearson's correlation coefficient and p-value
xtxt = min(DATA(:,1))+0.75*(max(DATA(:,1))-min(DATA(:,1)));
ytxt = 0.5*min(DATA(:,2))+0.5*max(DATA(:,2));
[R,P] = corrcoef(DATA(:,1),DATA(:,2))
if P(2,1)<0.001
    txt = sprintf('r = %0.2f*', R(2,1));
else
    txt = sprintf('r = %0.2f', R(2,1));
end
text(xtxt,ytxt,txt)

end