function fval = CalBDiff(t,vec,nx,dx)

% Parameters
Db = 20; %micrometer^2/s, Diffusion constant, CalB
alpha = Db/dx^2;
dbdx0 = 0; %flux db/dx at x=0
dbdxL = 0; %flux db/dx at x=L
btot(1:nx,1) = 40*10^-15; %micromole/micrometer^3, total CalB concentration
kbNeg = 19; %sec^-1; CalB used 
kbPos = 27*10^15; %micrometer^3/micromol s; CalB free

C = vec(1:nx,1);
b = vec(nx+1:2*nx,1);
    
% Stencil with Neuman boundary condition 
stencil = stencilMakerCalcium(nx);
stencil = alpha*stencil;

%Solve for dCdt
dbdt = stencil*b+(kbNeg*(btot-b)-kbPos*b.*C);

%First/last entry has constant flux from/to LHS/RHS 
dbdt(1,1) = dbdt(1,1) + 2*dx*dbdx0*alpha;
dbdt(nx,1) = dbdt(nx,1) - 2*dx*dbdxL*alpha; 

%To ensure size consistency for ode solvers
dbdt(nx+1:2*nx,1) = dbdt;

dbdt(1:nx,1) = 0;

fval = dbdt;






