function [CcNew,bNew,CeNew,pNew,beNew,RyRvalues,jP,jN,jR,jI,jS,jle,jlp,jERM,jSOC,jPM] = coupledEquations(tNew,tOld,speciesValues,RyRvalues,stencilOriginal,R,r,meanEdge,ptLst,iv)
%function [CcNew,bNew,CeNew,pNew,RyRvalue] = coupledEquations(tNew,tOld,speciesValues,nx,RyRvalues)
Cc = speciesValues(:,1);
b = speciesValues(:,2);
Ce = speciesValues(:,3);
p = speciesValues(:,4);
be = speciesValues(:,5);

if tOld < 1000
    jPM(1:iv.nx,1) = 0;
    jERM(1:iv.nx,1) = 0;
    jP(1:iv.nx,1) = 0;
    jN(1:iv.nx,1) = 0;
    jS(1:iv.nx,1) = 0;
    jI(1:iv.nx,1) = 0;
    jR(1:iv.nx,1) = 0;
    jSOC(1:iv.nx,1) = 0;
    poI(1:iv.nx,1) = 0;
    jle(1:iv.nx,1) = 0;
    jlp(1:iv.nx,1) = 0;
else
    [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,iv.nx,Cc,Ce,RyRvalues,iv.rhoR,iv);
    [jS,jP,jN,jle,jlp,jSOC] = pumpsEquations(Cc,Ce,iv);
    [jI, poI] = IP3REquations(Cc,Ce,p,iv);
    RyRvalues = [o1 o2 c1 c2];
    jPM = - jP - jN + jlp;
    jERM = jI + jR - jS + jle;
    
end

vec = [Cc;b;Ce;p;be];

[~, ySol] = Euler("FE",iv.nx,@(t,vec) diffusionEquations(t,vec,jERM,jPM,jSOC,stencilOriginal,R,r,meanEdge,ptLst,iv), [tOld,tNew], vec, Cc, iv);

CcNew = ySol(end,1:iv.nx).';
bNew = ySol(end,iv.nx+1:2*iv.nx).';
CeNew = ySol(end,2*iv.nx+1:3*iv.nx).';
pNew = ySol(end,3*iv.nx+1:4*iv.nx).';
beNew = ySol(end,4*iv.nx+1:5*iv.nx).';
end



