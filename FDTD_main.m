clear all
clc

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c_light = 1/sqrt(eps0*mu0);
freq = 2.0e9;

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
