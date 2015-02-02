% Construct a domain in three dimensions
n = 64;
domain = 0.2*[-n/2:n/2-1];
[X, Y, Z] = meshgrid(domain,domain,domain);

% Construct 3D gaussian centered at the origin.
gx = exp(-(domain-0).^2);
gy = exp(-(domain-0).^2);
gz = exp(-(domain-0).^2);
[GX, GY, GZ] = meshgrid(gx, gy, gz);
gaussian3D = GX.*GY.*GZ;

% Use isosurface to visualize the 3D gaussian. Each step visualizes the
% part of the gaussian at a certain value.
start = domain(1);
stop = domain(end);
Nsteps = 10;
for j = 0:Nsteps
    close all;
    isovalue = j/Nsteps;
    isosurface(X,Y,Z,gaussian3D,isovalue); 
    axis([start stop start stop start stop]);
    title(['3D Gaussian, Isosurface Value = ' mat2str(isovalue)]);
    grid on;
    draw now;
end