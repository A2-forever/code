mean_mass=sum(mass,2)/3;
mean_Gd=sum(Gd,2)/3;
mean_Cd=sum(Cd,2)/3;

subplot(2,1,1)
mean_Gd1=smooth(mean_Gd);
semilogx(mean_mass,mean_Gd1)

subplot(2,1,2)
semilogx(mean_mass,mean_Cd)