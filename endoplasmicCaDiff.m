function fval = endoplasmicCaDiff(t,Ce,nx,dx,jERM)

% Parameters
R = 0.4*10^-5;%0.2; %micrometer, Radius of Dendrite
r = 0.15*10^-5;%40*10^-3; %micrometer, Radius of ER

Dc = 220*10^-14;%220; %micrometer^2/s, Diffusion constant, IP3
alpha = Dc/dx^2;
dCdx0 = 0; %flux dp/dx at x=0
dCdxL = 0; %flux dp/dx at x=L 

% Stencil with Neuman boundary condition 
stencil = stencilMakerCalcium(nx);
stencil = alpha*stencil;

%Solve for dCdt
dCdt = stencil*Ce - 2*jERM/r;

%First/last entry has constant flux from/to LHS/RHS 
dCdt(1,1) = dCdt(1,1) + 2*dx*dCdx0*alpha;
dCdt(nx,1) = dCdt(nx,1) - 2*dx*dCdxL*alpha; 

fval = dCdt;






