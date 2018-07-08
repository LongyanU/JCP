clear;
clc;
close all

load('SaltA.mat')
aa=seis_record;
figure;imagesc(aa,[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')

load('SaltB.mat')
bb=seis_record;
figure;imagesc(bb,[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')

figure;imagesc(aa-bb,[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
