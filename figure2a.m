clear;
clc
close all

% [ 1.46951, -0.117903] 已经取得了很不错的结果

% 程序解释 p(1)+p(-1)-2p(0) =2cos(kh)-2   =   A
% 程序解释 p(2)+p(-2)-2p(0) =2cos(2kh)-2  =B

% coeff(1)*A + coeff(2)*B= -k^2*h^2;
% coeff是需要求解的系数



M=7;

AA=zeros(M,M);
b=zeros(M,1);
dt=0.002
dt=0.001

v=1500;
h=20;
k=linspace(1/50,0.8*pi/h,M);

for ii=1:M
    for kk=1:5
        xita=(kk-1)*pi/16;
        for jj=1:M
            %             AA(ii,jj)=2*cos(jj*k(ii)*h*cos(xita))-2 + 2*cos(jj*k(ii)*h*sin(xita))-2 +AA(ii,jj);
            AA(ii,jj)=2*cos(jj*k(ii)*h*cos(xita))-2*cos((jj-1)*k(ii)*h*cos(xita))...
                +2*cos(jj*k(ii)*h*sin(xita))-2*cos((jj-1)*k(ii)*h*sin(xita)) +AA(ii,jj);
        end
        b(ii)=1/(v^2*dt^2/h^2)*(2*cos(v*k(ii)*dt)-2)+b(ii);
    end
end

digits(6)
c=AA\b
vpa(c)'
tau=dt;

k=linspace(1/10000,pi/h,100);
r=v*dt/h;
temp=0;
for kk=1:5
    xita=(kk-1)*pi/16;
    temp=0;
    for m=1:M
        temp=temp+c(m)*(cos(m*k*h*cos(xita))-cos((m-1)*k*h*cos(xita))...
            +cos(m*k*h*sin(xita))-cos((m-1)*k*h*sin(xita)));
    end

    temp=1+temp*r^2;
    
    temp=acos(temp)./(k*v*tau);
    a1=(h/v*(1./temp-1));
    if kk==1
        figure;plot(a1,'k','LineWidth',2)
    elseif kk==2
        hold on;plot(a1,'m--','LineWidth',2);
    elseif kk==3
        hold on;plot(a1,'r:','LineWidth',2)
    elseif kk==4
        hold on; plot(a1,'b-.','LineWidth',2)
    else
        hold on;plot(a1,'c:.','LineWidth',2)
    end
    
end
legend('θ=0', 'θ=π/16','θ=2π/16','θ=3π/16','θ=4π/16')
xlabel('percentage of kh')
ylabel('\epsilon (\theta)')
grid on
axis([0 100 -1.5*10^-5  5*10^-5])