clear;
clc;
close all
%pseudo spectrum method May 10
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
izs=60;
ixs=nx/2;

%define source wavelet
f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff(diff(src))/dx^2;				% time derivative to obtain gaussian

seis_record=zeros(nt,nx);


%Initial wavefield
unow=zeros(nz,nx);             %Initializing current wavefield
uthen=zeros(nz,nx);            %Initialzing previous wavefield
cresult=zeros(nz,nx);          %Initialzing the migration result
factor=(dt.*velpart).^2;


for i=1:nx
    kx(i)=(2*pi/(nx*dx-dx)*(i-nx/2-1));
end

for i=1:nz
    kz(i)=(2*pi/(nz*dz-dz)*(i-nz/2-1));
end


%wavenumber calculation
for i=1:nz
    for j=1:nx
        spec(i,j)=power(kx(j),2)+power(kz(i),2);
    end
end

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

tic
for it=1:nt-2
    
    specxz=fft2(unow);
    %Shift to the center
    specxzshift=fftshift(specxz);
    
    %wavenumber shift
    specd2ud=specxzshift.*spec;
    sum = ifft2(ifftshift(specd2ud));
    
    
    %Computation of wavefield for current time step
    utemp=2*unow-factor.*(sum)-uthen;
    unow(izs,ixs)=unow(izs,ixs)+src(it)*dt^2;
    
    uthen=unow;
    unow=utemp;
    result(it,:)=unow(izs,:);
    
    a(it)=unow(nz/2-13,nx/2);
    b(it)=unow(nz/2+13,nx/2);
        [unow,uthen]=spongeABC(unow,uthen,nx,nz,45,45,0.009);
    
%     unow=taper.*unow;
%    uthen=taper.*uthen;
    
    if rem(it,isnap)== 0,
        imagesc(unow)
        axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(unow))))
        drawnow
    end
    
    
    
end
toc
figure;plot(-(a),'linewidth',2);hold on;plot(-(b),'r','linewidth',2)
xlabel('time(ms)')
ylabel('Amp')
legend('receiver A','receiver B')
grid on

save figure5c.mat
