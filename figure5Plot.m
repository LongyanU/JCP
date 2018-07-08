clear;
 close all;
 clc

load('figure5a.mat')
coeff
figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on
axis([ 0 1500 -7*10^-5 10*10^-5])

load('figure5b.mat')
figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
coeff
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on
axis([ 0 1500 -7*10^-5 10*10^-5])

load('figure5c.mat')
coeff
figure;plot(-(a),'linewidth',2);hold on;plot(-(b),'r','linewidth',2)
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on
axis([ 0 1500 -7*10^-5 10*10^-5])

load('figure5d.mat')
coeff
figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on
axis([ 0 1500 -7*10^-5 10*10^-5])

% load('figure5e.mat')
% coeff
% figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
% xlabel('time(ms)')
% ylabel('Amp')
% legend('receiver A','receiver B')
% grid on
% axis([ 0 1500 -7*10^-5 10*10^-5])
