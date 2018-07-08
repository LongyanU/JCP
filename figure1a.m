clear;clc;
close all
M=7;
h=20;
v=1500;
tau=0.002;
tau=0.0015;
r=v*tau/h;
% r=0;
%%%%%%%%%%%%%%%%%%%初始系数
c=zeros(1,M);
for m=1:M
    c(m)=1;
    temp=(-1)^(m+1)/(2*m-1);
    for n=1:M
        if n~=m
            c(m)=abs(((2*n-1)^2-r^2)/((2*n-1)^2-(2*m-1)^2))*c(m);
        end
    end
    c(m)=temp*c(m);
end
%%%%%%%%%%%%%%%%%%%%开始循环
r=v*tau/h;
k=linspace((pi)/(130*h),(pi)/h,120);
% k=linspace((pi)/(110*h),(pi)/h,110);
for loop=1:7
    AA=zeros(M,M);
    for i=1:100
        for kk=1:M
            a=0;atemp11=0;
            for jj=1:4
                atemp=0;
                xita=(jj-1)*15*pi/180;
                for n=1:M
                    atemp=atemp+c(n)*sin((n-0.5)*k(i)*h*cos(xita));
                end
                a=a+atemp^2;
                atemp11=atemp11+atemp*2*sin((kk-0.5)*k(i)*h*cos(xita));
            end
            
            b=0;btemp11=0;
            for jj=1:4
                btemp=0;
                xita=(jj-1)*15*pi/180;
                for n=1:M
                    btemp=btemp+c(n)*sin((n-0.5)*k(i)*h*sin(xita));
                end
                b=b+btemp^2;
                btemp11=btemp11+btemp*2*sin((kk-0.5)*k(i)*h*sin(xita));
            end
            
            g=1-1/2*r^2*(a+b);
            temp=1/(k(i)*tau)*1/sqrt(1-g^2)*1/2*r^2*(atemp11+btemp11);
            
            for j=1:M
                a=0;atemp11=0;
                for jj=1:4
                    atemp=0;
                    xita=(jj-1)*15*pi/180;
                    for n=1:M
                        atemp=atemp+c(n)*sin((n-0.5)*k(i)*h*cos(xita));
                    end
                    a=a+atemp^2;
                    atemp11=atemp11+atemp*2*sin((j-0.5)*k(i)*h*cos(xita));
                end
                
                b=0;btemp11=0;
                for jj=1:4
                    btemp=0;
                    xita=(jj-1)*15*pi/180;
                    for n=1:M
                        btemp=btemp+c(n)*sin((n-0.5)*k(i)*h*sin(xita));
                    end
                    b=b+btemp^2;
                    btemp11=btemp11+btemp*2*sin((j-0.5)*k(i)*h*sin(xita));
                end
                
                g=1-1/2*r^2*(a+b);
                tempvv=1/(k(i)*tau)*1/sqrt(1-g^2)*1/2*r^2*(atemp11+btemp11);
                
                AA(kk,j)=temp*tempvv+AA(kk,j);
            end
            
        end
    end
    
    ee=zeros(M,1);
    for i=1:100
        for kk=1:M
            
            a=0;atemp11=0;
            for jj=1:4
                atemp=0;
                xita=(jj-1)*15*pi/180;
                for n=1:M
                    atemp=atemp+c(n)*sin((n-0.5)*k(i)*h*cos(xita));
                end
                a=a+atemp^2;
                atemp11=atemp11+atemp*2*sin((kk-0.5)*k(i)*h*cos(xita));
            end
            
            b=0;btemp11=0;
            for jj=1:4
                btemp=0;
                xita=(jj-1)*15*pi/180;
                for n=1:M
                    btemp=btemp+c(n)*sin((n-0.5)*k(i)*h*sin(xita));
                end
                b=b+btemp^2;
                btemp11=btemp11+btemp*2*sin((kk-0.5)*k(i)*h*sin(xita));
            end
            
            g=1-1/2*r^2*(a+b);
            temp=1/(k(i)*tau)*1/sqrt(1-g^2)*1/2*r^2*(atemp11+btemp11);
            
            vv=acos(g)/(k(i)*tau);
            ee(kk)=temp*(vv-v)+ee(kk);
            
        end
    end
    
    x=-(AA\ee)
    
    for i=1:M
        c(i)=c(i)+x(i);
    end
    
end

%显示二维交错网格的网格频散
r=v*tau/h;

k=linspace((pi)/(100*h),(pi)/h,100);
for kk=1:5
    %     xita=(kkk-1)*15*pi/180;
    % xita=0
    % xita=pi/8
    %     xita=pi/4;
    xita=(kk-1)*pi/16;
    a=0;
    for n=1:M
        a=a+c(n)*sin((n-0.5)*k*h*cos(xita));
    end
    
    b=0;
    for n=1:M
        b=b+c(n)*sin((n-0.5)*k*h*sin(xita));
    end
    
    a=a.^2;
    b=b.^2;
    
    delta=2./(r*k*h);
    delta=delta.*asin(r*sqrt(a+b));
    a1=(h/v*(1./delta-1));
    %     a1=a1*10^5;
    %     hold on
    %     plot(v*k*h/(2*pi*h),(a1),'k')
    
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
grid on
axis([0 100 -1.5*10^-5  5*10^-5])
legend('θ=0', 'θ=π/16','θ=2π/16','θ=3π/16','θ=4π/16')
xlabel('percentage of kh')
ylabel('\epsilon (\theta)')
grid on