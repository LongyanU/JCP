clear;
clc;
close all
load('SaltA.mat')
figure;plot( seis_recordVx (:,400)+6*10^-7,'b','linewidth',2);


load('SaltB.mat')
hold on;plot( seis_recordVx (:,400)-6*10^-7,'r','linewidth',2);




% hold on;plot( aa-bb-5*10^-7,'m','linewidth',2);
% hold on;plot( dd-aa+5*10^-7,'m','linewidth',2);

xlabel('time(ms)')
ylabel('Amp')
title('')

% axis([0 4000 -0.82*10^-6 1.4*10^-6])
% hold on;plot(-seis_record(1500,:),'k')
yxis([-10^-3 10^-3])

legend('tra FD, LS coeff','Simplified FD, Linear coeff')
grid on
colormap gray