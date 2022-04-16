classdef getInitialValues
    % Obtain inital values/Equilibrium values 
    % for RyR probabilites and leakage velocities
    
    properties
        % Constants
        nx = 1;
        dt = 1;
        
        % Ca and CalB dynamics constants
        kbNeg = 19*10^-4; %10^-4 sec^-1; CalB used 
        kbPos = 27*10^-4; %micrometer^3/micromol s; CalB free
        
        % Initial/Equilibrium concentrations constants
        Ceq = 0.05; %mM
        b = 40; %micromole/micrometer^3, total CalB concentration, btot;
        Ce = 250; %micromole/micrometer^3
        p = 0*40*10^-3; %micromole/micrometer^3
        pr = 0*40*10^-3;
        be = 1000;
        
        % Pumps constants 
        Co = 1000; 
        CeEq = 250;
        
        % SERCA constants
        IS = 16.5*10^-19;%6.5*10^-21*10^-15; %mol micromol micrometer^-3 s^-1
        KS = 180*10^-3;%180*10^-18; %micromol/micrometer^3
        rhoS = 2390*10^10;%2390*10^10; %micrometer^-2
        
        % PMCA constants
        IP = 1.7*10^-21;%1.7*10^-23; %mol s^-1
        KP = 60*10^-3;%60*10^-18; %micromol/micrometer^3
        rhoP = 500*10^10;%500; %microm^-2
        
        % NCX constants
        IN = 2.5*10^-19;%2.5*10^-21; %mol s^-1
        KN = 1.8;%1.8*10^-15; %micromol micrometer^-3
        rhoN = 15*10^10;%15; %microm^-2
        
        %SOC constants
        vsoc = 10^-7;
        IORef = 2.1*10^-4*10^-15;
        z = 2;
        Ao = 2.5*10^-17;
        F = 96485;
        
        %Leakage Parameters
        vle = 3.7852160867136625448634251232620483162971680712871602736413478851318359375e-11;
        vlp = 4.497344938136755317628391160245617170533594109116393156000413000583648681640625e-12;

        % Diffusion constants
        Dc = 220*10^-14;%220; %micrometer^2/s, Diffusion constant
        Db = 20*10^-14;%20; %micrometer^2/s, Diffusion constant, CalB
        Dp = 280*10^-14;%280; %micrometer^2/s, Diffusion constant, IP3
        kp = 0.11*10^-1;%0.11; %sec^-1;  
        
        % RyR Probabilities constants
        prob = [0 0 0 0];
        kapos = 1500*10^-4;
        kaneg = 28.8*10^-4;
        kbpos = 1500*10^-4;
        kbneg = 385.9*10^-4;
        kcpos = 1.75*10^-4;
        kcneg = 0.1*10^-4;
        IRRef = 3.5*10^-16;%3.5*10^-18*10^6*10^-15; %micromole s^-1
        rhoR = 3*10^10;%3.0; %micrometer^-2
        CeRef = 250;%250*10^-15; %micromole micrometer^-3
        
        % IP3REquations constants
        IIRef = 1.1*10^-17;%1.1*10^-19*10^6; %micromol s^-1
        rhoI = 173*10^9;%17.3; %micrometer^-2
        d1 = 0.13;%0.13*10^-15; %micromol/microm^3
        d2 = 1.05;%1.05*10^-15; %micromol/microm^3
        d3 = 0.94;%0.94*10^-15; %micromol/microm^3
        d5 = 82.3*10^-3;%82.3*10^-18; %micromol/microm^3

    end
    methods
        function RyRProb = getInitProb(obj,Cc)
            % This returns intital/equilibirum probabilites for RyR ODEs
            % for <Cc> concentration
            A = [1 1 1 1; obj.kaneg 0 -obj.kapos*Cc^4 0; obj.kbpos*Cc^3 -obj.kbneg 0 0; obj.kcpos 0 0 -obj.kcneg];
            B = [10^6; 0; 0; 0];
            RyRProb = round(linsolve(A,B));
        end
        function [vle,vlp] = getvlevlp(obj,Cc,Ce,Co,p,RyRvalues)
            tOld = 1;
            tNew = 10;
            [jR,o1,o2,c1,c2] = RyREquations(tOld,tNew,1,Cc,Ce,RyRvalues,obj.rhoR,obj);
            [jS,jP,jN,jle,jlp,jSOC] = pumpsEquations(Cc,Ce,obj);
            [jI, poI] = IP3REquations(Cc,Ce,p,obj);
            vlp = (jP + jN)/(Co - Cc);
            vle = (- jI - jR + jS)/(Ce - Cc);
        end
    end 
end