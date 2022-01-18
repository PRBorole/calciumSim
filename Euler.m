function [t,vecSol] = Euler(type,ddt,tSpan,dt,vec,Cc)

if type =="FE"  
    tStart = tSpan(1);
    tEnd = tSpan(2);
    t = tStart:dt:tEnd-dt;
    j=1;
    % vecSol = zeros(size(t,2),size(vec,1));
    vecSol = zeros(1,size(vec,1));
    for i = t
        fn = ddt(i,vec);
        vec = vec + fn*dt;
        vecSol(1,:) = vec;
        % vecSol(j,:) = vec;
        j=j+1;
    end
else
    % RyR Probabilities constants
    nx = length(Cc);
    kaPos = 1500*10^-4;
    kaNeg = 28.8*10^-4;
    kbPos = 1500*10^-4;
    kbNeg = 385.9*10^-4;
    kcPos = 1.75*10^-4;
    kcNeg = 0.1*10^-4;
    tStart = tSpan(1);
    tEnd = tSpan(2);
    t = tStart:dt:tEnd;
    % vecSol = zeros(size(t,2),size(vec,1));
    vecSol = zeros(1,size(vec,1));
    
    for i = t
        i = (i-tStart)/dt+1;
        for k = 1:nx
            A = [-(kaNeg+kbPos*Cc(k)^3+kcPos) kbNeg kaPos*Cc(k)^4 kcNeg; kbPos*Cc(k)^3 -kbNeg 0 0; kaNeg 0 -kaPos*Cc(k)^4 0; kcPos 0 0 -kcNeg];
            v = [vec(k); vec(k+nx); vec(k+2*nx); vec(k+3*nx)];
            v = linsolve((eye(4)-A*dt),v);
            vecSol(1,k) = v(1);
            vecSol(1,k+nx) = v(2);
            vecSol(1,k+2*nx) = v(3);
            % vecSol(i,k+3*nx) = v(4);
            vecSol(1,k+3*nx) = 10^6 - v(1) - v(2) - v(3);
        end
        vec = vecSol(1,:);
    end
end
end
