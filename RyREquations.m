%function [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,Cc,Ce,RyRvalues,rhoR)
function [jR,o1,o2,c1,c2,tSol,vecSol] = RyREquations(tOld,tNew,Cc,Ce,RyRvalues,rhoR)

IRRef = 3.5*10^-16;%3.5*10^-18*10^6*10^-15; %micromole s^-1
CeRef = 250;%250*10^-15; %micromole micrometer^-3
%rhoR = 3*10^10;%3.0; %micrometer^-2
nx = length(Cc);

o1 = RyRvalues(:,1);
o2 = RyRvalues(:,2);
c1 = RyRvalues(:,3);
c2 = RyRvalues(:,4);

vec = [o1;o2;c1;c2];

%options = odeset('RelTol',1e-2,'AbsTol',1e-1);

% [tSol, vecSol] = ode45(@(t,vec) do1c1c2o2dt(t,vec,Cc), [tOld,tNew], vec);
% dt = 1;
% [tSol, vecSol] = Euler("FE",@(t,vec) do1c1c2o2dt(t,vec,Cc), [tOld,tNew],dt,vec,Cc);
dt = 1;
[tSol, vecSol] = Euler("BE",@(t,vec) do1c1c2o2dt(t,vec,Cc), [tOld,tNew],dt,vec,Cc);
% [tSol, vecSol] = ode15s(@(t,vec) do1c1c2o2dt(t,vec,Cc), [tOld,tNew], vec);

vecSol = round(vecSol);
o1 = vecSol(end,1:nx).';
o2 = vecSol(end,nx+1:2*nx).';
c1 = vecSol(end,2*nx+1:3*nx).';
c2 = vecSol(end,3*nx+1:4*nx).';

poR = o1+o2;

IR = IRRef*(Ce-Cc)/CeRef;
jR = rhoR*poR.*IR;
jR = jR./10^6;
end

function fval = do1c1c2o2dt(t,vec,Cc)

nx = length(Cc);
o1 = vec(1:nx,1);
o2 = vec(nx+1:2*nx,1);
c1 = vec(2*nx+1:3*nx,1);
c2 = 10^6 -o1 -o2 -c1;
% c2 = vec(3*nx+1:4*nx,1);


% dc1dt
kaPos = 1500*10^-4;%1500 * 10^60; %micrometer^12 micromole^-4 s^-1
kaNeg = 28.8*10^-4;%28.8; %s^-1
dc1dt = kaNeg*o1 - kaPos*(Cc.^4).*c1;

% do2dt
kbPos = 1500*10^-4;%1500 * 10^45; %micrometer^9 micromol^-3 s^-1
kbNeg = 385.9*10^-4;%385.9; %s^-1
do2dt = kbPos*(Cc.^3).*o1 - kbNeg*o2;

% dc2dt
kcPos = 1.75*10^-4;%1.75; %s^-1
kcNeg = 0.1*10^-4;%0.1; %s^-1
dc2dt = kcPos*o1 - kcNeg*c2;

% do1dt
do1dt = -dc1dt -dc2dt -do2dt;
fval = [do1dt;do2dt;dc1dt;dc2dt];
end