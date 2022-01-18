function fval = diffusionEquations(t,vec,nx,jERM,jPM,jSOC,stencilOriginal,R,r,meanEdge,ptLst)

Dc = 220*10^-14;%220; %micrometer^2/s, Diffusion constant
Db = 20*10^-14;%20; %micrometer^2/s, Diffusion constant, CalB

btot(1:nx,1) = 40;%40*10^-15; %micromole/micrometer^3, total CalB concentration
kbNeg = 19*10^-4;%19; %sec^-1; CalB used 
kbPos = 27*10^-4;%27*10^15; %micrometer^3/micromol s; CalB free

C = vec(1:nx,1);
b = vec(nx+1:2*nx,1);
Ce = vec(2*nx+1:3*nx,1);
p = vec(3*nx+1:4*nx,1);

nTrig = length(ptLst);
[jSYNPCa,jSYNPp] = inputSynapse(t,stencilOriginal,R,r,meanEdge,ptLst);

% Cytosol Ca2+
% Stencil with Neuman boundary condition 
stencil = stencilOriginal;
stencil = Dc*stencil;

%Solve for dCdt
dCdt = stencil*C+kbNeg*(btot-b)-kbPos*b.*C+2*r.*jERM./(R.^2 - r.^2)+2*R.*jPM./(R.^2 - r.^2);
for i = 1:nTrig
    dCdt(ptLst(i),1) = dCdt(ptLst(i),1) + jSYNPCa(i);
end
    
% CalBindin
% Stencil with Neuman boundary condition 
stencil = stencilOriginal;
stencil = Db*stencil;

%Solve for dbdt
dbdt = stencil*b+kbNeg*(btot-b)-kbPos*b.*C;

% ER Ca2+
% Stencil with Neuman boundary condition 
stencil = stencilOriginal;
stencil = Dc*stencil;

%Solve for dCdt
dCedt = stencil*Ce - 2*jERM./r + 2*jSOC./r;

% IP3
Dp = 280*10^-14;%280; %micrometer^2/s, Diffusion constant, IP3
pr(1:nx,1) = 40*10^-3;%40*10^-18; %micromol/micrometer^3, Basal IP3 concentration
kp = 0.11*10^-1;%0.11; %sec^-1;  

% Stencil with Neuman boundary condition 
stencil = stencilOriginal;
stencil = Dp*stencil;

%Solve for dpdt
dpdt = stencil*p-kp*(p-pr);
for i = 1:nTrig
    dpdt(ptLst(i),1) = dpdt(ptLst(i),1) + jSYNPp(i);
end

fval = [dCdt; dbdt; dCedt; dpdt];
end

function [jSYNPCa,jSYNPp] = inputSynapse(t,stencilOriginal,R,r,meanEdge,ptLst)

nTrig = length(ptLst);
tStartCa = 1500;
tEndCa = 1510;
tStartp = 1500;
tEndp = 3500;

%jIN1 = 2.849925377*10^-6;

jINCa = 0*2.5*10^-6;
slopeCa = -jINCa/10;

jINp = 5*10^-6;
slopep = -jINp/2000;

jSYNPCa(1:nTrig) = 0;
jSYNPp(1:nTrig) = 0;

for i = 1:nTrig
    pt = ptLst(i);
    flag = isEnd(pt,stencilOriginal);
    if t>tStartCa && t<tEndCa
        jINCay = slopeCa*(t-tStartCa)+jINCa;
        if flag
            jSYNPCa(i) = 2*jINCay/meanEdge;
        else
            jSYNPCa(i) = 2*R(pt)*jINCay/(R(pt)^2-r(pt)^2);
        end
    else
        jSYNPCa(i) = 0;
    end
    
    if t>tStartp && t<tEndp
        jINpy = slopep*(t-tStartp)+jINp;
        if flag
            jSYNPp(i) = 2*jINpy/meanEdge;
        else
            jSYNPp(i) = 2*R(pt)*jINpy/(R(pt)^2-r(pt)^2);
        end
    else
        jSYNPp(i) = 0;
    end
end
end

function flag = isEnd(pt,stencilOriginal)
% Determines if selected point is end point
if nnz(stencilOriginal(pt,:)) < 3
    flag = 1;
else
    flag = 0;
end
end 
