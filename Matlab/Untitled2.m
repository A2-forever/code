logcell=log(cell);
p1=polyfit(logcell(2:end,1),logcell(2:end,2),1);
y1=polyval(p1,log(cell(2:end,1)))';
hold on
plot(logcell(2:end,1),logcell(2:end,2),'r')
plot(logcell(2:end,1),y1,'b')
legend('计算结果','拟合曲线:y=0.6942x+1.4208')
xlabel('logTIME')
ylabel('logMSD')
hold off