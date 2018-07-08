%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unsplit complex frequency shifted perfectly matched layer for second-order wave equation using auxiliary differential equations
% Gaoet al.: JASA Express Letters
%%%%%%%%%                  simple model                     %%%%%%%%%%%%%%%%%
clear;
clc;
close all
%混合偏导数用伪谱法，空间偏导数用有限差分 May 10
% Elapsed time is 3.984338 seconds.
nt=1520;
dx=20;
dz=20;
nx=250;
nz=250;


v=ones(nz,nx)*2500;
v(1:nz/2,:)=1500;
isnap=60;    % snapshot sampling


%define the axes
x=(0:nx-1)*dx;
z=(0:nz-1)*dz;
dt=0.0015;

velpart=v;


%define source wavelet
f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff(diff(src))/dx^2;				% time derivative to obtain gaussian

seis_record=zeros(nt,nx);

temp=[ 1.57866, -0.296598, 0.0949307, -0.0344762, 0.0120067, -0.00344529, 0.000605554];

coeff=[-temp(1)*2  temp(1)-temp(2) temp(2)-temp(3) temp(3)-temp(4) temp(4)-temp(5) temp(5)-temp(6) temp(6)-temp(7) temp(7)];


%%%%%%%%%%%%%%%%%%%%pml%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=12*dx;
xita=pi/8;
pd=2;
Vmax=max(max(v));
R=10^-3;
d0=(pd+1)*Vmax*log(1/R)/(2*delta);

ddx=zeros(nz,nx);
ddz=zeros(nz,nx);
%left side
for i=1:nz
    for j=20:-1:9
        temp=dx*(20-j);
        ddx(i,j)=d0*(temp/delta)^2;
    end
end

%right side
for i=1:nz
    for j=nx-20:nx-9
        temp=dx*(j-nx+20);
        ddx(i,j)=d0*(temp/delta)^2;
    end
end

%top side
for i=20:-1:9
    for j=1:nx
        temp=dx*(20-i  +1);
        ddz(i,j)=d0*(temp/delta)^2;
    end
end

%bottom side
for i=nz-20:nz-9
    for j=1:nx
        temp=dx*(i-nz+20  +1);
        ddz(i,j)=d0*(temp/delta)^2;
    end
end

ddx2=ddx.*ddx;
ddz2=ddz.*ddz;
deriDx=zeros(nz,nx);
deriDz=zeros(nz,nx);
%left side
for i=1:nz
    for j=20:-1:9
        temp=dx*(20-j);
        deriDx(i,j)=d0*2*temp*(1/delta)^2;
    end
end
%right side
for i=1:nz
    for j=nx-20:nx-9
        temp=dx*(j-nx+20);
        deriDx(i,j)=d0*2*temp*(1/delta)^2;
    end
end
%top side
for i=20:-1:9
    for j=1:nx
        temp=dx*(20-i  +1);
        deriDz(i,j)=d0*2*temp*(1/delta)^2;
    end
end

%bottom side
for i=nx-20:nx-9
    for j=1:nx
        temp=dx*(i-nz+20  +1);
        deriDz(i,j)=d0*2*temp*(1/delta)^2;
    end
end
alfax=zeros(nz,nx);
alfaz=zeros(nz,nx);
alfaMax=f0;
%left side
for i=1:nz
    for j=20:-1:9
        temp=dx*(20-j);
        alfax(i,j)=alfaMax*(1-temp/delta);
    end
end
%right side
for i=1:nz
    for j=nx-20:nx-9
        temp=dx*(j-nx+20);
        alfax(i,j)=alfaMax*(1-temp/delta);
    end
end
%top side
for i=20:-1:9
    for j=1:nx
        temp=dx*(20-i +1);
        alfaz(i,j)=alfaMax*(1-temp/delta);
    end
end

%bottom side
for i=nz-20:nz-9
    for j=1:nx
        temp=dx*(i-nz+20 +1);
        alfaz(i,j)=alfaMax*(1-temp/delta);
    end
end

deriAlfa=alfaMax*(-1/delta);
%%%%%%%%%%%%%%%%%%%%pml%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seis_record=zeros(nt,nx);
tic
u3=zeros(nz,nx);
u2=zeros(nz,nx);
u1=zeros(nz,nx);

%Z direction
u3z=zeros(nz,nx);
u2z=zeros(nz,nx);
u1z=zeros(nz,nx);
dz=dx;
d1px=zeros(nz,nx);
d1pz=zeros(nz,nx);

p=zeros(nz,nx);
d2px=p; d2pz=p;
pnew=p;pold=p;
addpoints=p;
% Source location
zs=nx/2;
xs=nz/2;

zs=60;
xs=nx/2;

M=7;
for it=1:nt-2,
    for j=9:nz-9,
        for k=9:nx-9,
            
            if k <=20
                d1px(j,k)=(p(j,k)-p(j,k+1)); % x direction
            else
                d1px(j,k)=(p(j,k)-p(j,k-1)); % x direction
            end
            if j<=20
                d1pz(j,k)=(p(j,k)-p(j+1,k)); % z direction
            else
                d1pz(j,k)=(p(j,k)-p(j-1,k)); % z direction
            end
            
            d2pz(j,k)=coeff(1)*p(j,k);
            for ii=2:M+1
                d2pz(j,k)=d2pz(j,k)+coeff(ii)*(p(j-ii+1,k)+p(j+ii-1,k));
            end
            
            d2px(j,k)=coeff(1)*p(j,k);
            for ii=2:M+1
                d2px(j,k)=d2px(j,k)+coeff(ii)*(p(j,k-ii+1)+p(j,k+ii-1));
            end
            
%             addpoints(j,k)=2*coeff(M+2)*( p(j+1,k+1)+p(j+1,k-1)+p(j-1,k+1)+p(j-1,k-1)-p(j+1,k)-p(j-1,k)-p(j,k+1)-p(j,k-1) );
        end
    end
    %x direction
    u3=u3+dt*(-ddx.*u3-alfax.*u3+ddx2.*(deriDx+deriAlfa).*d1px/dx);
    u2=u2+dt*(-ddx.*u2-alfax.*u2+ ddx.*(2*deriDx.*d1px/dx+deriAlfa.*d1px/dx+ddx.*d2px/(dx^2)) - u3);
    u1=u1+dt*(-ddx.*u1-alfax.*u1+deriDx.*d1px/dx+2*ddx.*d2px/(dx^2)-u2);
    
    %z direction
    u3z=u3z+dt*(-ddz.*u3z-alfaz.*u3z+ddz2.*(deriDz+deriAlfa).*d1pz/dz);
    u2z=u2z+dt*(-ddz.*u2z-alfaz.*u2z+ ddz.*(2*deriDz.*d1pz/dz+deriAlfa.*d1pz/dz+ddz.*d2pz/(dz^2)) - u3z);
    u1z=u1z+dt*(-ddz.*u1z-alfaz.*u1z+deriDz.*d1pz/dz+2*ddz.*d2pz/(dx^2)-u2z);
    
    pnew=2*p-pold+v.*v.*(d2px+d2pz+addpoints)*dt^2/dx^2 -u1.*v.*v.*dt^2  -u1z.*v.*v.*dt^2;
    pnew(zs,xs)=pnew(zs,xs)+src(it)*dt^2;
    pold=p;											% time lev(k,j)els
    p=pnew;
    
    a(it)=p(nz/2-13,nx/2);
    b(it)=p(nz/2+13,nx/2);
    
    if rem(it,isnap)== 0,
        imagesc(p)
        axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' time step: %i - max ampl: %g ',it,max(max(p))))
        drawnow
    end
  
    
    if it==600
        pp1=p;
    elseif it==720
        pp2=p;
    elseif it==740
        pp3=p;
    elseif it==780
        pp4=p;
    end
end
toc
figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on
axis([ 0 1500 -7*10^-5 10*10^-5])
save figure5d.mat