function fval = IP3Diff(t,p,nx,dx)

% Parameters
R = 0.4*10^-5;%0.2; %micrometer, Radius of Dendrite
r = 0.15*10^-5;%40*10^-3; %micrometer, Radius of ER

Dp = 280*10^-14;%280; %micrometer^2/s, Diffusion constant, IP3
alpha = Dp/dx^2;
dpdx0 = 0; %flux dp/dx at x=0
dpdxL = 0; %flux dp/dx at x=L
pr(1:nx,1) = 40*10^-3;%40*10^-18; %micromol/micrometer^3, Basal IP3 concentration
kp = 0.11*10^-4;%0.11; %sec^-1;  
 
% Stencil with Neuman boundary condition 
stencil = stencilMakerCalcium(nx);
stencil = alpha*stencil;

%Solve for dpdt
dpdt = stencil*p-kp*(p-pr)./(R.^2 - r.^2);

%First/last entry has constant flux from/to LHS/RHS 
dpdt(1,1) = dpdt(1,1) + 2*dx*dpdx0*alpha;
dpdt(nx,1) = dpdt(nx,1) - 2*dx*dpdxL*alpha; 

fval = dpdt;






