% Elapsed time is 496.169502 seconds.seconds.
% May 16

% Elapsed time is 764.150298 seconds.  June 29
% 时间已过 781.736860 秒。 秒。 July 7,2018.
clear
clc %%%%%%%
close all
nt=4020;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')

c1=flipud(c);

load('Figure7aCoeff.mat')
v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end


clear v
v=vv;
% v=ones(nz,nx)*1500;
% v(1:nz/2,:)=1500;

coeff=zeros(nz,nx,7);
for ii=1:nz
    for jj=1:nx
        temp111=floor((v(ii,jj)-1486)/10)+1;
        coeff(ii,jj,:)=coeffJune(temp111 ,:);
    end
end


dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=600-150+25;

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;
p=zeros([nz nx]);

r=v*dt/h;

x11=d2px;  z11=d2pz;

taper=ones(nz,nx);
for i=1:50
    for j=1:nx
        taper(i,j)=0.5-0.5*cos(pi*(i-1)/(50-1));
        taper(nz-i+1,j)=taper(i,j);
    end
end
for i=1:nz
    for j=1:50
        taper(i,j)=taper(i,j)*(0.5-0.5*cos(pi*(j-1)/(50-1)));
        taper(i,nx-j+1)=taper(i,j);
    end
end



p=zeros([nz nx]); Vx=p; Vz=p;

seis_recordVz=zeros(nt,nx);
seis_recordVx=zeros(nt,nx);
tic
for it=1:nt-2,
    
    d2pz11=Vz-circshift(Vz,[0 1]);
    d2pz12=(circshift(Vz,[0 -1])-circshift(Vz,[0 2]));
    d2pz13=(circshift(Vz,[0 -2])-circshift(Vz,[0 3]));
    d2pz14=(circshift(Vz,[0 -3])-circshift(Vz,[0 4]));
    d2pz15=(circshift(Vz,[0 -4])-circshift(Vz,[0 5]));
    d2pz16=(circshift(Vz,[0 -5])-circshift(Vz,[0 6]));
    d2pz17=(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    
    d2px11=Vx-circshift(Vx,[1 0]);
    d2px12=(circshift(Vx,[-1 0])-circshift(Vx,[2 0]));
    d2px13=(circshift(Vx,[-2 0])-circshift(Vx,[3 0]));
    d2px14=(circshift(Vx,[-3 0])-circshift(Vx,[4 0]));
    d2px15=(circshift(Vx,[-4 0])-circshift(Vx,[5 0]));
    d2px16=(circshift(Vx,[-5 0])-circshift(Vx,[6 0]));
    d2px17=(circshift(Vx,[-6 0])-circshift(Vx,[7 0]));
    
    
    %     d2px=coeff(1)*d2px11+coeff(2)*d2px12+coeff(3)*d2px13+coeff(4)*d2px14+coeff(5)*d2px15+coeff(6)*d2px16...
    %         +coeff(7)*d2px17;
    %     d2pz=coeff(1)*d2pz11+coeff(2)*d2pz12+coeff(3)*d2pz13+coeff(4)*d2pz14+coeff(5)*d2pz15+coeff(6)*d2pz16...
    %         +coeff(7)*d2pz17;
    
    
    d2px=coeff(:,:,1).*d2px11+coeff(:,:,2).*d2px12+coeff(:,:,3).*d2px13+coeff(:,:,4).*d2px14+coeff(:,:,5).*d2px15+coeff(:,:,6).*d2px16...
        +coeff(:,:,7).*d2px17;
    d2pz=coeff(:,:,1).*d2pz11+coeff(:,:,2).*d2pz12+coeff(:,:,3).*d2pz13+coeff(:,:,4).*d2pz14+coeff(:,:,5).*d2pz15+coeff(:,:,6).*d2pz16...
        +coeff(:,:,7).*d2pz17;
    
    p=p-dt*v.^2.*(d2px+d2pz)/h;
    
    p(zs,xs)= p(zs,xs)+src(it)*dt^2;
    seis_record(it,:)=p(zs,:);
    
    d2pz1=(circshift(p,[0 -1])-circshift(p,[0 0]));
    d2pz2=(circshift(p,[0 -2])-circshift(p,[0 1]));
    d2pz3=(circshift(p,[0 -3])-circshift(p,[0 2]));
    d2pz4=(circshift(p,[0 -4])-circshift(p,[0 3]));
    d2pz5=(circshift(p,[0 -5])-circshift(p,[0 4]));
    d2pz6=(circshift(p,[0 -6])-circshift(p,[0 5]));
    d2pz7=(circshift(p,[0 -7])-circshift(p,[0 6]));
    
    
    
    d2px1=(circshift(p,[-1])-circshift(p,[0]));
    d2px2=(circshift(p,[-2])-circshift(p,[1]));
    d2px3=(circshift(p,[-3])-circshift(p,[2]));
    d2px4=(circshift(p,[-4])-circshift(p,[3]));
    d2px5=(circshift(p,[-5])-circshift(p,[4]));
    d2px6=(circshift(p,[-6])-circshift(p,[5]));
    d2px7=(circshift(p,[-7])-circshift(p,[6]));
    %     d2px=coeff(1)*d2px1+coeff(2)*d2px2+coeff(3)*d2px3+coeff(4)*d2px4+coeff(5)*d2px5+coeff(6)*d2px6...
    %         +coeff(7)*d2px7;
    %     d2pz=coeff(1)*d2pz1+coeff(2)*d2pz2+coeff(3)*d2pz3+coeff(4)*d2pz4+coeff(5)*d2pz5+coeff(6)*d2pz6...
    %         +coeff(7)*d2pz7;
    d2px=coeff(:,:,1).*d2px1+coeff(:,:,2).*d2px2+coeff(:,:,3).*d2px3+coeff(:,:,4).*d2px4+coeff(:,:,5).*d2px5+coeff(:,:,6).*d2px6...
        +coeff(:,:,7).*d2px7;
    d2pz=coeff(:,:,1).*d2pz1+coeff(:,:,2).*d2pz2+coeff(:,:,3).*d2pz3+coeff(:,:,4).*d2pz4+coeff(:,:,5).*d2pz5+coeff(:,:,6).*d2pz6...
        +coeff(:,:,7).*d2pz7;
    
    Vx=Vx-dt*d2px/h;
    Vz=Vz-dt*d2pz/h;
    
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
        [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
%     Vx=taper.*Vx;
%     Vz=Vz.*taper; % for reviewing
    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    
    
    
    if it==500
        pp1=p;
    elseif it==1000
        pp2=p;
    elseif it==1500
        pp3=p;
    elseif it==2000
        pp4=p;
    elseif it==2500
        pp5=p;
    end
    
end
toc

figure;imagesc(v)
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs,zs,'*r')

save('SaltA.mat')
