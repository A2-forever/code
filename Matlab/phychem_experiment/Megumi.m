x=((elec(1,:)-elec))./(T);

p1=polyfit(x(2:end-1,1),elec(2:end-1,1),1);
y1=polyval(p1,x(2:end-1,1))';

p2=polyfit(x(2:end,2),elec(2:end,2),1);
y2=polyval(p2,x(2:end,2))';

% hold on
% plot(x(2:end-1,1),elec(2:end-1,1),"b*");
% plot(x(2:end-1,1),y1,"*b-");
% hold off


hold on
plot(x(2:end,2),elec(2:end,2),"r*");
plot(x(2:end,2),y2,"r*-");
hold off