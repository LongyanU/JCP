clear;
clc
close all

M=7;

AA=zeros(M,M);
b=zeros(M,1);


v=1500;
h=15;
for loop=1:3
    if loop==1
        M=3;
        k=linspace(1/50,0.35*pi/h,M);
        iii=1;
        for dt=0.00005:0.00005:0.01;
            AA=zeros(M,M);
            b=zeros(M,1);
            rr(iii)=v*dt/h;
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
            
            
            c=AA\b;
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
        figure; plot(rr,(rr),'r','linewidth',2)
        hold on;plot(rr,(ss),'k','linewidth',2);
    elseif loop==2
        M=5;
        k=linspace(1/50,0.55*pi/h,M);
        
        iii=1;
        for dt=0.00005:0.00005:0.01;
            AA=zeros(M,M);
            b=zeros(M,1);
            rr(iii)=v*dt/h;
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
            
            
            c=AA\b;
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
        
        hold on;plot(rr,(ss),'b','linewidth',2);
    else
        M=7;
        k=linspace(1/50,0.78*pi/h,M);
        iii=1;
        for dt=0.00005:0.00005:0.01;
            AA=zeros(M,M);
            b=zeros(M,1);
            rr(iii)=v*dt/h;
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
            
            
            c=AA\b;
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
         hold on;plot(rr,(ss),'m','linewidth',2);
    end
    
   
end
legend('r','M=3','M=5','M=7')
xlabel('r')
ylabel('Stability')
grid on
