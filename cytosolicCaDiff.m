function fval = cytosolicCaDiff(t,vec,nx,dx)

% Parameters
Dc = 220; %micrometer^2/s, Diffusion constant
alpha = Dc/dx^2;

if t < 1*10^-3 %For pulse of 1 ms
    dCdx0 = 2.5*10^-12; %flux dC/dx at x=0, micromole micrometer^-2 s^-1
else
    dCdx0 = 0;
end

dCdxL = 0; %flux dC/dx at x=L
btot(1:nx,1) = 40*10^-15; %micromole/micrometer^3, total CalB concentration
kbNeg = 19; %sec^-1; CalB used 
kbPos = 27*10^15; %micrometer^3/micromol s; CalB free

C = vec(1:nx,1);
b = vec(nx+1:2*nx,1);

% Stencil with Neuman boundary condition 
stencil = stencilMakerCalcium(nx);
stencil = alpha*stencil;

%Solve for dCdt

dCdt = stencil*C+(kbNeg*(btot-b)-kbPos*b.*C);

%First/last entry has constant flux from/to LHS/RHS 
dCdt(1,1) = dCdt(1,1) - 2*dx*dCdx0*alpha;
dCdt(nx,1) = dCdt(nx,1) + 2*dx*dCdxL*alpha; 

%To ensure size consistency for ode solvers
dCdt(nx+1:2*nx,1) = 0;

fval = dCdt;






