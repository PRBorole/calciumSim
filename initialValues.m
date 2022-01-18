function fval = initialValues(Cc,nx)


for density=[0:5*10^9:5*10^10]
    density
    [jR,o1,o2,c1,c2] = RyREquations(0,10^6,load("Cc200.mat").Ceq,250,[324 0 994014 10^6-994014],density);
    jI = IP3REquations(load("Cc200.mat").Ceq,250,0.04);
    [jS,jP,jN,jle,jlp] = pumpsEquations(load("Cc200.mat").Ceq,250);
    num2str((jS-jR-jI)/(250-load("Cc200.mat").Ceq),100)
end

% kaPos = 1500*10^-3;
% kaNeg = 28.8*10^-3;
% kbPos = 1500*10^-3;
% kbNeg = 385.9*10^-3;
% kcPos = 1.75*10^-3;
% kcNeg = 0.1*10^-3;
% 
% A = [1 1 1 1; kaNeg 0 -kaPos*Cc^4 0; kbPos*Cc^3 -kbNeg 0 0; kcPos 0 0 -kcNeg];
% B = [1; 0; 0; 0];
% RyRProb = linsolve(A,B);
% o1(1:nx,1) = RyRProb(1);
% o2(1:nx,1) = RyRProb(2);
% c1(1:nx,1) = RyRProb(3);
% c2(1:nx,1) = 1-o1-o2-c1;
% fval = [o1 o2 c1 c2];
% 
% p = 0.04; % microM
% Ce = 250; %microM
% Co = 2000;
% 
% %SERCA Parameters
% IS = 6.5*10^-18;%6.5*10^-21*10^-15; %mol micromol micrometer^-3 s^-1
% KS = 180*10^-3;%180*10^-18; %micromol/micrometer^3
% rhoS = 2390*10^10;%2390; %micrometer^-2
% jS = rhoS*IS*Cc./((KS+Cc).*Ce);
% 
% %PMCA Parameters
% IP = 1.7*10^-20;%1.7*10^-23; %mol s^-1
% KP = 60*10^-3;%60*10^-18; %micromol/micrometer^3
% rhoP = 500*10^10;%500; %microm^-2
% jP = rhoP*IP*Cc.^2./(KP^2+Cc.^2);
% 
% %NCX Parameters
% IN = 2.5*10^-18;%2.5*10^-21; %mol s^-1
% KN = 1.8;%1.8*10^-15; %micromol micrometer^-3
% rhoN = 15*10^10;%15; %microm^-2
% jN = rhoN*IN*Cc./(KN+Cc);
% 
% vlp = (jP + jN)/(Co - Cc);
% save("vlp.mat","vlp");
% 
% [jR,o1,o2,c1,c2] = RyREquations(0,100,Cc,Ce,[RyRProb(1) RyRProb(2) RyRProb(3) RyRProb(4)]);
% jI = IP3REquations(Cc,Ce,p);
% 
% vle = (-jI - jR + jS)/(Ce - Cc);
% save("vle.mat","vle");

end
