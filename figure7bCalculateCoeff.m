clear;
clc
close all

load ('vv.mat') % then we obtain veloicty model c

M=7;

AA=zeros(M,M);
b=zeros(M,1);
dt=0.001;


v=1500;
h=20;
k=linspace(1/50,0.8*pi/h,M);

coeffJune28=zeros(floor(max(max(c))- min(min(c)) )+1, M);
tic;
iii=1;

for v= floor(min(min(c))):floor(max(max(c)))
    for ii=1:M
        for kk=1:5
            xita=(kk-1)*pi/16;
            for jj=1:M
                AA(ii,jj)=2*cos(jj*k(ii)*h*cos(xita))-2*cos((jj-1)*k(ii)*h*cos(xita))...
                    +2*cos(jj*k(ii)*h*sin(xita))-2*cos((jj-1)*k(ii)*h*sin(xita)) +AA(ii,jj);
            end
            b(ii)=1/(v^2*dt^2/h^2)*(2*cos(v*k(ii)*dt)-2)+b(ii);
        end
    end
    
    coeffJune28(iii,:)=AA\b;
    iii=iii+1;
end
toc
save('Figure7bCoeff.mat')