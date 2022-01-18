function [jS,jP,jN,jle,jlp,jSOC] = pumpsEquations(Cc,Ce)

Co = 1000;%2; %mM = 2*10^3*10^-15; %to micromol/micrometer^3
CeEq = 250;

%SERCA Parameters
IS = 6.5*10^-19;%6.5*10^-21*10^-15; %mol micromol micrometer^-3 s^-1
KS = 180*10^-3;%180*10^-18; %micromol/micrometer^3
rhoS = 2390*10^10;%2390*10^10; %micrometer^-2
jS = rhoS*IS*Cc./((KS+Cc).*Ce);

%PMCA Parameters
IP = 1.7*10^-21;%1.7*10^-23; %mol s^-1
KP = 60*10^-3;%60*10^-18; %micromol/micrometer^3
rhoP = 500*10^10;%500; %microm^-2
jP = rhoP*IP*Cc.^2./(KP^2+Cc.^2);

%NCX Parameters
IN = 2.5*10^-19;%2.5*10^-21; %mol s^-1
KN = 1.8;%1.8*10^-15; %micromol micrometer^-3
rhoN = 15*10^10;%15; %microm^-2
jN = rhoN*IN*Cc./(KN+Cc);

%SOC
vsoc = 10^-10;
jSOC = vsoc.*(CeEq - Ce);


%Leakage Parameters

vle = 3.7852160867136625448634251232620483162971680712871602736413478851318359375e-11;
vlp = 4.497344938136755317628391160245617170533594109116393156000413000583648681640625e-12;

jle = vle.*(Ce-Cc);
jlp = vlp.*(Co-Cc);

end