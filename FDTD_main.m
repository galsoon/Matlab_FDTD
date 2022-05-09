clear all
close all
clc

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c_light = 1/sqrt(eps0*mu0);
freq = 2.0e9;
len_tf = 0.5;
number_of_materials = 2;

L=[2.5 2.5];
nx = 500;
ny = 500;
X=linspace(0,L(1),nx)-L(1)/2;
Y=linspace(0,L(2),ny)-L(2)/2;

dx = X(nx) - X(nx-1);
dy = Y(ny) - Y(ny-1);
dt = (1/c_light/sqrt(1/(dx^2)+1/(dy^2)))*0.99;

Index = zeros(nx,ny);
IndexX = zeros(nx,ny);
IndexY = zeros(nx,ny);

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
k_Ez_b = 2.0/(2.0*eps0*Material(1,1)+Material(1,3)*dt);
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

k_Fz_1 = zeros(ny,1);
k_Fz_2 = zeros(ny,1);
k_Ez_1 = zeros(nx,1);
k_Ez_2 = zeros(nx,1);
k_Gx_1 = zeros(ny-1,1);
k_Gx_2 = zeros(ny-1,1);
k_Hx_1 = zeros(nx,1);
k_Hx_2 = zeros(nx,1);
k_Gy_1 = zeros(nx-1,1);
k_Gy_2 = zeros(nx-1,1);
k_Hy_1 = zeros(ny,1);
k_Hy_2 = zeros(ny,1);
k_Hz_1 = zeros(nx,1);
k_Hz_2 = zeros(nx,1);
k_Ex_1 = zeros(nx,1);
k_Ex_2 = zeros(nx,1);
k_Ey_1 = zeros(ny,1);
k_Ey_2 = zeros(ny,1);

k_Fz_1(:) = 1.0;
k_Fz_2(:) = dt;
k_Ez_1(:) = 1.0;
k_Ez_2(:) = 1.0/eps0;
k_Gx_1(:) = 1.0;
k_Gx_2(:) = dt/dy;
k_Hx_1(:) = 1.0/mu0;
k_Hx_2(:) = 1.0/mu0;
k_Gy_1(:) = 1.0;
k_Gy_2(:) = dt/dx;
k_Hy_1(:) = 1.0/mu0;
k_Hy_2(:) = 1.0/mu0;
k_Hz_1(:) = 1.0;
k_Hz_2(:) = 1.0/mu0;
k_Ex_1(:) = 2.0;
k_Ex_2(:) = 2.0;
k_Ey_1(:) = 2.0;
k_Ey_2(:) = 2.0;

K_a1(1:number_of_materials) = 1./(2*eps0*Material(1:number_of_materials,1)+Material(1:number_of_materials,3)*dt);
K_b1(1:number_of_materials) = 2*eps0*Material(1:number_of_materials,1)-Material(1:number_of_materials,3)*dt;

k_Fz_2 = repmat(k_Fz_2(2:ny-1)',nx-2,1);
k_Gx_1 = repmat(k_Gx_1(1:ny-1)',nx,1);
k_Gx_2 = repmat(k_Gx_2(1:ny-1)',nx,1);
k_Gy_1 = repmat(k_Gy_1(1:nx-1),1,ny);
k_Gy_2 = repmat(k_Gy_2(1:nx-1),1,ny);
k_Hz_2 = repmat(k_Hz_2(2:nx-1),1,ny-2);
k_Ey_1 = repmat(k_Ey_1(1:ny)',nx-1,1);
k_Ey_2 = repmat(k_Ey_2(1:ny)',nx-1,1);
k_Ez_1 = repmat(k_Ez_1(2:nx-1),1,ny-2);
k_Ez_2 = repmat(k_Ez_2(2:nx-1),1,ny-2);


figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'doublebuffer','on');


for T=1:200
%    TE Mode ******************************************************************************************* 
%     Hz의 초기 생성 위치

%     Ey_1D(1:nx_b) = k_Ey_a(1:nx_b).*Ey_1D(1:nx_b) - k_Ey_b(1:nx_b).*(Hz_1D(2:nx_b+1)-Hz_1D(1:nx_b));
    Ey_1D(1:399) = k_Ey_a(1:399).*Ey_1D(1:399) - k_Ey_b(1:399).*(Hz_1D(2:400)-Hz_1D(1:399));
    Hz_1D(1) = sin(2*pi*freq*(T-1)*dt);
%     Hz_1D(2:nx_b) = Hz_1D(2:nx_b) - k_Hz_b(2:nx_b).*(Ey_1D(2:nx_b)-Ey_1D(1:nx_b-1));
    Hz_1D(2:399) = Hz_1D(2:399) - k_Hz_b(2:399).*(Ey_1D(2:399)-Ey_1D(1:398));

    
%     Hz 처리

%     Fz_r1 = Wz(2:nx-1,2:ny-1);
    Fz_r1 = Wz(2:499,2:499);
%     Wz(2:nx-1,2:ny-1) = Wz(2:nx-1,2:ny-1) + k_Fz_2.*((Ex(2:nx-1,2:ny-1)-Ex(2:nx-1,1:ny-2))/dy - (Ey(2:nx-1,2:ny-1)-Ey(1:nx-2,2:ny-1))/dx);
    Wz(2:499,2:499) = Wz(2:499,2:499) + k_Fz_2.*((Ex(2:499,2:499)-Ex(2:499,1:498))/dy - (Ey(2:499,2:499)-Ey(1:498,2:499))/dx);
%     Hz(2:nx-1,2:ny-1) = Hz(2:nx-1,2:ny-1) + k_Hz_2.*(Wz(2:nx-1,2:ny-1)-Fz_r1);
    Hz(2:499,2:499) = Hz(2:499,2:499) + k_Hz_2.*(Wz(2:499,2:499)-Fz_r1);
    
    
%     Hz 위줄, 아랫줄에서 반사되는 wave 제거

%     Hz(nx_a+1,ny_a+1:ny_b+1) = Hz(nx_a+1,ny_a+1:ny_b+1) + dt./(mu0*Material(Index(nx_a+1,ny_a+1:ny_b+1)+1,2)'*dx)*Ey_1D(1);
    Hz(101,101:400) = Hz(101,101:400) + dt./(mu0*Material(Index(101,101:400)+1,2)'*dx)*Ey_1D(1);
%     Hz(nx_b+1,ny_a+1:ny_b+1) = Hz(nx_b+1,ny_a+1:ny_b+1) - dt./(mu0*Material(Index(nx_b+1,ny_a+1:ny_b+1)+1,2)'*dx)*Ey_1D(nx_b-nx_a+2);
    Hz(400,101:400) = Hz(400,101:400) - dt./(mu0*Material(Index(400,101:400)+1,2)'*dx)*Ey_1D(301);

    
%     Ex 처리

%     Gx_r1 = Mx(1:nx,1:ny-1);
    Gx_r1 = Mx(1:500,1:499);
%     Mx(1:nx,1:ny-1) = k_Gx_1.*Mx(1:nx,1:ny-1) + k_Gx_2.*(Hz(1:nx,2:ny)-Hz(1:nx,1:ny-1));
    Mx(1:500,1:499) = k_Gx_1.*Mx(1:500,1:499) + k_Gx_2.*(Hz(1:500,2:500)-Hz(1:500,1:499));
%     Ex(1:nx,1:ny-1) = K_a1(IndexX(1:nx,1:ny-1)+1).*(K_b1(IndexX(1:nx,1:ny-1)+1).*Ex(1:nx,1:ny-1)+k_Ex_1.*Mx(1:nx,1:ny-1)-k_Ex_2.*Gx_r1);
    Ex(1:500,1:499) = K_a1(Index(1:500,1:499)+1).*(K_b1(IndexX(1:500,1:499)+1).*Ex(1:500,1:499)+k_Ex_1.*Mx(1:500,1:499)-k_Ex_2.*Gx_r1);
    
    
%     Ey 처리

%     Gy_r1 = My(1:nx-1,1:ny);
    Gy_r1 = My(1:499,1:500);
%     My(1:nx-1,1:ny) = k_Gy_1.*My(1:nx-1,1:ny) - k_Gy_2.*(Hz(2:nx,1:ny)-Hz(1:nx-1,1:ny));
    My(1:499,1:500) = k_Gy_1.*My(1:499,1:500) - k_Gy_2.*(Hz(2:500,1:500)-Hz(1:499,1:500));
%     Ey(1:nx-1,1:ny) = K_a1(IndexY(1:nx-1,1:ny)+1).*(K_b1(IndexY(1:nx-1,1:ny)+1).*Ey(1:nx-1,1:ny)+k_Ey_1.*My(1:nx-1,1:ny)-k_Ey_2.*Gy_r1);
    Ey(1:499,1:500) = K_a1(IndexY(1:499,1:500)+1).*(K_b1(IndexY(1:499,1:500)+1).*Ey(1:499,1:500)+k_Ey_1.*My(1:499,1:500)-k_Ey_2.*Gy_r1);
    
    
%     좌,우측에서 반사되는 wave 처리 

%     Ex(nx_a+1:nx_b+1,ny_a) = Ex(nx_a+1:nx_b+1,ny_a) - 2*dt/dy*K_a1(Index(nx_a+1:nx_b+1,ny_a)+1)'.*Hz_1D(2:(nx_b-nx_a+2));
    Ex(101:400,100) = Ex(101:400,100) - 2*dt/dy*K_a1(Index(101:400,100)+1)'.*Hz_1D(2:301);
%     Ex(nx_a+1:nx_b+1,ny_b+1) = Ex(nx_a+1:nx_b+1,ny_b+1) + 2*dt/dy*K_a1(Index(nx_a+1:nx_b+1,ny_b+1)+1)'.*Hz_1D(2:(nx_b-nx_a+2));
    Ex(101:400,400) = Ex(101:400,400) + 2*dt/dy*K_a1(Index(101:400,400)+1)'.*Hz_1D(2:301);
    
    
%     위,아래에서 반사되는 wave처리

%     Ey(nx_a,ny_a+1:ny_b+1) = Ey(nx_a,ny_a+1:ny_b+1) + 2*dt/dx*K_a1(Index(nx_a,ny_a+1:ny_b+1)+1,1)'.*Hz_1D(2);
    Ey(100,101:400) = Ey(100,101:400) + 2*dt/dx*K_a1(Index(100,101:400)+1,1)'.*Hz_1D(2);
%     Ey(nx_b+1,ny_a+1:ny_b+1) = Ey(nx_b+1,ny_a+1:ny_b+1) - 2*dt/dx*K_a1(Index(nx_b+1,ny_a+1:ny_b+1)+1,1)'.*Hz_1D(nx_b-nx_a+2);
    Ey(400,101:400) = Ey(400,101:400) - 2*dt/dx*K_a1(Index(400,101:400)+1,1)'.*Hz_1D(301);
    
    
    
    
 

    
    if(mod(T,10) == 0)
        
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
        drawnow;
    end
end

