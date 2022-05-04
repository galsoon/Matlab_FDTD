clear all
clc

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c_light = 1/sqrt(eps0*mu0);
freq = 2.0e9;
len_tf = 0.5;

L=[2.5 2.5];
nx = 500;
ny = 500;
X=linspace(0,L(1),nx)-L(1)/2;
Y=linspace(0,L(2),ny)-L(2)/2;

dx = X(nx) - X(nx-1);
dy = Y(ny) - Y(ny-1);
dt = (1/c_light/sqrt(1/(dx^2)+1/(dy^2)))*0.99;

Index = zeros(nx,ny);
Index_X = zeros(nx,ny);
Index_Y = zeros(nx,ny);

coout = 100;
Material = zeros(2,3);
Material(1,1) = 1.0;
Material(1,2) = 1.0;
Material(1,3) = 0;
Material(2,1) = 1.0;
Material(2,2) = 1.0;
Material(2,3) = 1e70;

nx_a = round(len_tf/dx);
nx_b = round((L(1)-len_tf)/dx);
ny_a = round(len_tf/dy);
ny_b = round((L(2)-len_tf)/dy);

Ez_1D = zeros(nx_b+1,1);
Fz_1D = zeros(nx_b+1,1);
Hy_1D = zeros(nx_b,1);
k_Fz_a = zeros(nx_b+1,1);
k_Fz_b = zeros(nx_b+1,1);
k_Hy_a = zeros(nx_b,1);
k_Hy_b = zeros(nx_b,1);
k_Ez_a = (2.0*eps0*Material(1,1)-Material(1,3)*dt)/(2.0*eps0*Material(1,1)+Material(1,3)*dt);
Hz_1D = zeros(nx_b+1,1);
Ey_1D = zeros(nx_b,1);
k_Hz_a = zeros(nx_b+1,1);
k_Hz_b = zeros(nx_b+1,1);
k_Ey_a = zeros(nx_b,1);
k_Ey_b = zeros(nx_b,1);

k_Hy_a(:) = 1.0;
k_Hy_b(:) = dt/(mu0*Material(1,2)*dx);
k_Fz_a(:) = 1.0;
k_Fz_b(:) = dt/dx;
k_Hz_a(:) = 1.0;
k_Hz_b(:) = dt/(mu0*Material(1,2)*dx);
k_Ey_a(:) = 1.0;
k_Ey_b(:) = dt/(eps0*Material(1,1)*dx);

k_a_max = 1;
m = 4;
R_err = 1e-16;
eta = sqrt(mu0*Material(1,2)/eps0/Material(1,1));

Fz = zeros(nx,ny);
Tz = zeros(nx,ny);
Gx = zeros(nx,ny-1);
Gy = zeros(nx-1,ny);
Ez = zeros(nx,ny);
Hx = zeros(nx,ny-1);
Hy = zeros(nx-1,ny);
Wz = zeros(nx,ny);
Mx = zeros(nx,ny-1);
My = zeros(nx-1,ny);
Hz = zeros(nx,ny);
Ex = zeros(nx,ny-1);
Ey = zeros(nx-1,ny);

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'doublebuffer','on');
tic


for i=1:5
    subplot(1,2,1);
    pcolor(Y(1:ny),X(1:nx),Ez(1:nx,1:ny));
    shading interp;
    axis image;
    colorbar;
    
    
    subplot(1,2,2);
    pcolor(Y(1:ny),X(1:nx),Hz(1:nx,1:ny));
    shading interp;
    axis image;
    colorbar;
end

