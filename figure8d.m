clear;
clc;
close all


load('SaltA.mat')
% plotimage(-seis_record)
% xlabel('x/dx')
% ylabel('time(ms)')
% title('')
hold on;plot(-seis_record(1:end,400),'b','linewidth',2)

aa=-seis_record(1:end,400);


load('SaltB.mat')
% plotimage(-seis_record)
% xlabel('x/dx')
% ylabel('time(ms)')
% title('')
hold on;plot(-seis_record(1:end,400),'k','linewidth',2)
bb=-seis_record(1:end,400);

cc=aa-bb;
hold on;plot(cc-1.5*10^-3,'c','linewidth',2)
% hold on;plot(real(seis_record(1500,:)),'m')
legend('tra FD, LS coeff','New FD, Linear coeff','the difference between the two')

axis([500 4000 -2.2*10^-3 3.7*10^-3])
xlabel('Travel time(ms)')
ylabel('Amp')
grid on
box on