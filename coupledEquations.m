function [CcNew,bNew,CeNew,pNew,RyRvalues,jP,jN,jR,jI,jS,poI] = coupledEquations(tNew,tOld,speciesValues,nx,RyRvalues,stencilOriginal,R,r,meanEdge,ptLst)
%function [CcNew,bNew,CeNew,pNew,RyRvalue] = coupledEquations(tNew,tOld,speciesValues,nx,RyRvalues)
Cc = speciesValues(:,1);
b = speciesValues(:,2);
Ce = speciesValues(:,3);
p = speciesValues(:,4);

rhoR = 3*10^10;
if tOld < 1000
    jPM(1:nx,1) = 0;
    jERM(1:nx,1) = 0;
    jP(1:nx,1) = 0;
    jN(1:nx,1) = 0;
    jS(1:nx,1) = 0;
    jI(1:nx,1) = 0;
    jR(1:nx,1) = 0;
    jSOC(1:nx,1) = 0;
    poI(1:nx,1) = 0;
else
    [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,Cc,Ce,RyRvalues,rhoR);
    [jS,jP,jN,jle,jlp,jSOC] = pumpsEquations(Cc,Ce);
    [jI, poI] = IP3REquations(Cc,Ce,p);
    RyRvalues = [o1 o2 c1 c2];
    jPM = - jP - jN + jlp;
    jERM = jI + jR - jS + jle;
end

vec = [Cc;b;Ce;p];
% options = odeset('RelTol',10^0,'AbsTol',10^0);
% [t, CcbSol] = ode15s(@(t,vec) diffusionEquations(t,vec,nx,jERM,jPM,stencilOriginal,R,r,meanEdge,ptLst), [tOld:1:tNew], vec,options);
% t
% (t(end)-t(1))/size(t,1)

dt = 1;
[~, ySol] = Euler("FE",@(t,vec) diffusionEquations(t,vec,nx,jERM,jPM,jSOC,stencilOriginal,R,r,meanEdge,ptLst), [tOld,tNew], dt, vec);

CcNew = ySol(end,1:nx).';
bNew = ySol(end,nx+1:2*nx).';
CeNew = ySol(end,2*nx+1:3*nx).';
pNew = ySol(end,3*nx+1:4*nx).';
end



