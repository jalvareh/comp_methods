% Initiate workspace.
clear all; close all; clc; 

% Load data.
load Testdata;

% Initialize objects.
L=15;        % spatial domain size
n=64;        % number of Fourier modes
T = 20;      % number of time points
totalk = zeros(n,n,n);

% Construct spatial domain.
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x; 
[X,Y,Z]=meshgrid(x,y,z); 

% Construct frequency domain.
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
mink = ks(1);
maxk = ks(end);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Sum the fourier transform of the 3D signal.
for j=1:T 
    Un(:,:,:) = reshape(Undata(j,:),n,n,n);
    totalk = totalk + fftn(Un);
%     close all, isosurface(X,Y,Z,abs(Un),0.4); 
%     axis([-20 20 -20 20 -20 20]), grid on, drawnow; 
%     pause(0.1); 
end 

% Calculate the average, normalized frequency spectrum.
avgk = abs((totalk))/T;

% Find the location of the maximum peak in frequency space, which is the 
% frequency of the marbel.
[maxK, index] = max(avgk(:));
[kyloc, kxloc, kzloc] = ind2sub(size(avgk),index);
normk = avgk/maxK;

% Use isosurface to visualize the main peak in the frequency spectrum.
% for j = 0:2:10
%     close all;
%     isosurface(Kx,Ky,Kz,normk,j/10); 
%     axis([mink maxk mink maxk mink maxk]);
%     title('FT');
%     grid on;
%     drawnow;
% end
% 
% % Construct 3D gaussian filter, centered about the peak.
w = 2;                                % Filter width.
fx = exp(-w*(ks-ks(kxloc)).^2);
fy = exp(-w*(ks-ks(kyloc)).^2);
fz = exp(-w*(ks-ks(kzloc)).^2);
[FX, FY, FZ] = meshgrid(fx, fy, fz);
filter = FX.*FY.*FZ;                  % Final 3D Gaussian filter.

% Use isosurface to visualize the filter in the frequency spectrum.
for j = 0:10
    close all;
    isosurface(Kx,Ky,Kz,fftshift(filter),j/10); 
    axis([mink maxk mink maxk mink maxk]);
    title('Filter');
    grid on;
    drawnow;
end
% 
% % Use isosurface to visualize the filter signal in the frequency spectrum.
% for j = 0:2:10
%     close all;
%     isosurface(Kx,Ky,Kz,normk.*filter,j/10); 
%     axis([mink maxk mink maxk mink maxk]);
%     title('Filtered FT');
%     grid on;
%     drawnow;
% end


% Loop through the data (over time), applying the filter to the 3D signal
% at each time.
% fprintf('x\ty\tz\n');
% close all;
xloc = zeros(1,T);
yloc = zeros(1,T);
zloc = zeros(1,T);
for j=1:T  
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    unk = fftn(Un);
    unkf = unk.*filter;
    x = abs(ifftn((unkf)));
    [dummy, index] = max(x(:));
    [yind, xind, zind] = ind2sub(size(x),index);
    xloc(j) = x2(xind);
    yloc(j) = x2(yind);
    zloc(j) = x2(zind);
end 

% close all;
figure(1);
subplot(3,1,1);
plot(xloc,'k.-', 'MarkerSize', 20);
xlabel('t'), ylabel('x');
title('Components of Marble Trajectory');
grid on;
subplot(3,1,2);
plot(yloc,'k.-', 'MarkerSize', 20);
xlabel('t'), ylabel('y');
grid on;
subplot(3,1,3);
plot(zloc,'k.-', 'MarkerSize', 20);
xlabel('t'), ylabel('z');
grid on;

figure(2);
plot3(xloc, yloc, zloc, 'b-', 'LineWidth', 2);
hold on;
plot3(xloc(end), yloc(end), zloc(end), 'rx', 'MarkerSize', 50,...
    'LineWidth',4);
grid on;
xlabel('x'), ylabel('y'), zlabel('z');
title('Marble Trajectory');

fprintf('Final location: x = %0.1f, y = %0.1f, z = %0.1f\n', ...
    xloc(end), yloc(end), zloc(end));

