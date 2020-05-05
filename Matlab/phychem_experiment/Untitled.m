ps=101.325;
p_d=data(:,1)/ps;
p_d(:,2)=p_d./((1-p_d).*data(:,2));
p=polyfit(p_d(1:end-2,1),p_d(1:end-2,2),1);
plot(p_d(1:end-2,1),p_d(1:end-2,2));