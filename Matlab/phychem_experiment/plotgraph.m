n=3;

hold on
for i=1:n
    plot(name{i}(:,1),name{i}(:,2)-(i-3)*1000)
end
legend("��Ʒa��XRDͼ","��Ʒb��XRDͼ","��Ʒc��XRDͼ")
hold off
% i=10;
% plot(name{i}(:,1),name{i}(:,2))