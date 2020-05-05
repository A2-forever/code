n=3;

hold on
for i=1:n
    plot(name{i}(:,1),name{i}(:,2)-(i-3)*1000)
end
legend("样品a的XRD图","样品b的XRD图","样品c的XRD图")
hold off
% i=10;
% plot(name{i}(:,1),name{i}(:,2))