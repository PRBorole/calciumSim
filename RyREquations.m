%function [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,Cc,Ce,RyRvalues,rhoR)
function [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,nx,Cc,Ce,RyRvalues,rhoR,iv)

o1 = RyRvalues(:,1);
o2 = RyRvalues(:,2);
c1 = RyRvalues(:,3);
c2 = RyRvalues(:,4);

vec = [o1;o2;c1;c2];

[tSol, vecSol] = Euler("BE",nx,@(t,vec) do1c1c2o2dt(t,vec,Cc), [tOld,tNew],vec,Cc,iv);

vecSol = round(vecSol);
o1 = vecSol(end,1:nx).';
o2 = vecSol(end,nx+1:2*nx).';
c1 = vecSol(end,2*nx+1:3*nx).';
c2 = vecSol(end,3*nx+1:4*nx).';

poR = o1+o2;
IR = iv.IRRef*(Ce-Cc)/iv.CeRef;
jR = iv.rhoR*poR.*IR;
jR = jR./10^6;
end

function fval = do1c1c2o2dt(t,vec,Cc)

nx = length(Cc);
o1 = vec(1:nx,1);
o2 = vec(nx+1:2*nx,1);
c1 = vec(2*nx+1:3*nx,1);
c2 = 10^6 -o1 -o2 -c1;

dc1dt = iv.kaneg*o1 - iv.kapos*(Cc.^4).*c1; % dc1dt
do2dt = iv.kbpos*(Cc.^3).*o1 - iv.kbneg*o2; % do2dt
dc2dt = iv.kcpos*o1 - iv.kcneg*c2; % dc2dt
do1dt = -dc1dt -dc2dt -do2dt; % do1dt
fval = [do1dt;do2dt;dc1dt;dc2dt];
end