function [jI, poI] = IP3REquations(Cc,Ce,p)
IIRef = 1.1*10^-17;%1.1*10^-19*10^6; %micromol s^-1
CeRef = 250;%250*10^-15; %micromole/micrometer^3
rhoI = 173*10^9;%17.3; %micrometer^-2
d1 = 0.13;%0.13*10^-15; %micromol/microm^3
d2 = 1.05;%1.05*10^-15; %micromol/microm^3
d3 = 0.94;%0.94*10^-15; %micromol/microm^3
d5 = 82.3*10^-3;%82.3*10^-18; %micromol/microm^3

poI = ((d2*Cc.*p)./((Cc.*p+d2*p+d3*Cc+d1*d2).*(Cc+d5))).^3;


II = IIRef*(Ce-Cc)/CeRef;
jI = rhoI*poI.*II;
end
