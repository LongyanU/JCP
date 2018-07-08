clear;
clc;
close all
load('SaltA.mat')
figure;imagesc(-seis_recordVx,[-10^-6 10^-6])
aa=-seis_recordVx;
xlabel('x/dx')
ylabel('time(ms)')
title('')
colormap gray

load('SaltB.mat')
figure;imagesc(-seis_recordVx,[-10^-6 10^-6])
bb=-seis_recordVx;
colormap gray
xlabel('x/dx')
ylabel('time(ms)')
title('')

cc=aa-bb;
figure;imagesc(cc,[-10^-6 10^-6])
bb=-seis_recordVx;
xlabel('x/dx')
ylabel('time(ms)')
title('')
colormap gray

