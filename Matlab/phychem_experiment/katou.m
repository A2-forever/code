subplot(2,1,1)
Gd1=ones(95,3);
for i=1:3
    Gd1(:,i)=smooth(Gd(:,i));
    semilogx(mass(:,i),Gd1(:,i))
    hold on
    for j=1:3
        if i==3 && j==3
            break
        end
%         text(S1(j,3*i-2),S1(j,3*i-1),S2(j,3*i-1))
    end
        
end
text
hold off

subplot(2,1,2)
for i=1:3
    semilogx(mass(:,i),Cd(:,i))
    hold on
    for j=1:3
        if i==3 && j==3
            break
        end
%         text(S1(j,3*i-2),S1(j,3*i),S2(j,3*i))
    end
end
hold off