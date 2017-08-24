% Test-case-1D Poisson Equation
clc; close all; clear all;
L = 1;											% Length of the domain
N = 10;											% Number of cells
for N = N:10:1000
x = (L/N)*([1:N] - 1/2);    % Cell centers
phi_an_1 = sin(2*pi*x);			% First analytical solution (F_AS)
rho = 4*pi^2*sin(2*pi*x);		% Source term corresponding to F_AS
phi = PoiSolv1D(rho,0,0,L/N);
% plot(x,phi_an_1,'r-',x,phi,'b-')
semilogy(N,max(abs(phi-phi_an_1)),'bo',N,(L/N)^2,'ro','LineWidth',2)
hold on
end
grid on