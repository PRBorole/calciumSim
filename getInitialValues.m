classdef getInitialValues
    % Obtain inital values/Equilibrium values 
    % for RyR probabilites and leakage velocities
    
    properties
        % RyR Probabilities constants
        kaPos = 1500*10^-4;
        kaNeg = 28.8*10^-4;
        kbPos = 1500*10^-4;
        kbNeg = 385.9*10^-4;
        kcPos = 1.75*10^-4;
        kcNeg = 0.1*10^-4;
        
        % Leakage velocity constants
    end
    methods
        function RyRProb = getInitProb(obj,Cc)
            % This returns intital/equilibirum probabilites for RyR ODEs
            % for <Cc> concentration
            A = [1 1 1 1; obj.kaNeg 0 -obj.kaPos*Cc^4 0; obj.kbPos*Cc^3 -obj.kbNeg 0 0; obj.kcPos 0 0 -obj.kcNeg];
            B = [1; 0; 0; 0];
            RyRProb = linsolve(A,B);
        end
        function [vle,vlp] = getvlevlp(obj,Cc,Ce,Co,p,RyRvalues,rhoR)
            tOld = 1;
            tNew = 10;
            [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,Cc,Ce,RyRvalues,rhoR);
            [jS,jP,jN,jle,jlp] = pumpsEquations(Cc,Ce);
            [jI, poI] = IP3REquations(Cc,Ce,p);
            vlp = (jP + jN)/(Co - Cc); 
            vle = (- jI - jR + jS)/(Ce - Cc);
        end
    end 
end