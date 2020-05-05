mean1=sum(Data,2);
mean1([1,5,6])=mean1([1,5,6])/3;
mean1([2,3,4])=mean1([2,3,4])/2;

yeta_ex=mean1(:)/mean1(1)*yeta;

X=[0,1,2/3,1/2,1/3,1/4]';
C_x=X*C;

yeta_sp=(yeta_ex-yeta_ex(1))/yeta_ex(1);
yeta_r=yeta_ex/yeta_ex(1);
        
yeta_sp_C=yeta_sp(2:end)./C_x(2:end);
lnyeta_r=log(yeta_r);
lnyeta_r_C=lnyeta_r(2:end)./C_x(2:end);

hold on
p1=polyfit(C_x(2:end),yeta_sp_C,1);
y1=polyval(p1,C_x)';
plot(C_x,y1)
scatter(C_x(2:end),yeta_sp_C,"x")

p2=polyfit(C_x(2:end),lnyeta_r_C,1);
y2=polyval(p2,C_x)';
plot(C_x,y2)
scatter(C_x(2:end),lnyeta_r_C,"x")
hold off

