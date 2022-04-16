function [jI, poI] = IP3REquations(Cc,Ce,p,iv)

poI = ((iv.d2*Cc.*p)./((Cc.*p+iv.d2*p+iv.d3*Cc+iv.d1*iv.d2).*(Cc+iv.d5))).^3;
II = iv.IIRef*(Ce-Cc)/iv.CeRef;
jI = iv.rhoI*poI.*II;
end
