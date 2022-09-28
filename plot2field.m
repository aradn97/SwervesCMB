close all

figure(1)
load -ascii clt.d
semilogx(clt(:,1)*(2.723e6)^2,'linewidth',1.5,'color',[100,121,162]/256)
hold on
semilogx(clt(:,2)*(2.723e6)^2,'linewidth',1.5,'color',[59,113,86]/256)
semilogx(clt(:,3)*(2.723e6)^2,'linewidth',1.5,'color',[135,81,81]/256)
semilogx(clt(:,4)*(2.723e6)^2,'linewidth',1.5,'color',[255,177,100]/256)
legend('(1/9)C_l^{adb}','(4/25)C_l^{iso}','(-4/15)C_l^{cross}','C_l^{total}');
ylabel('l(l+1)/2\pi C_l^{TT}')
xlabel('l')

figure(2)
load -ascii cle.d
semilogx(cle(:,1)*(2.723e6)^2,'linewidth',1.5,'color',[100,121,162]/256)
hold on
semilogx(cle(:,2)*(2.723e6)^2,'linewidth',1.5,'color',[59,113,86]/256)
semilogx(cle(:,3)*(2.723e6)^2,'linewidth',1.5,'color',[135,81,81]/256)
semilogx(cle(:,4)*(2.723e6)^2,'linewidth',1.5,'color',[255,177,100]/256)
legend('(1/9)C_l^{adb}','(4/25)C_l^{iso}','(-4/15)C_l^{cross}','C_l^{total}');
ylabel('l(l+1)/2\pi C_l^{EE}')
xlabel('l')

figure(3)
load -ascii clc.d
plot(clc(:,1),'linewidth',1.5,'b')
hold on
plot(clc(:,2),'linewidth',1.5,'r')
plot(clc(:,3),'linewidth',1.5,'g')
plot(clc(:,4),'linewidth',1.5,'k')
legend('(1/9)C_l^{adb}','(4/25)C_l^{iso}','(-4/15)C_l^{cross}','C_l^{total}');
ylabel('l(l+1)/2\pi C_l^{TE}')
xlabel('l')