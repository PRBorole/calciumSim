function [t,vecSol] = Euler(type,nx,ddt,tSpan,vec,Cc,iv)

if type =="FE"  
    tStart = tSpan(1);
    tEnd = tSpan(2);
    t = tStart:iv.dt:tEnd-iv.dt;
    j=1;
    % vecSol = zeros(size(t,2),size(vec,1));
    vecSol = zeros(1,size(vec,1));
    for i = t
        fn = ddt(i,vec);
        vec = vec + fn*iv.dt;
        vecSol(1,:) = vec;
        % vecSol(j,:) = vec;
        j=j+1;
    end
else
    % RyR Probabilities constants
    tStart = tSpan(1);
    tEnd = tSpan(2);
    t = tStart:iv.dt:tEnd;
    % vecSol = zeros(size(t,2),size(vec,1));
    vecSol = zeros(1,size(vec,1));
    
    for i = t
        i = (i-tStart)/iv.dt+1;
        for k = 1:nx
            A = [-(iv.kaneg+iv.kbpos*Cc(k)^3+iv.kcpos) iv.kbneg iv.kapos*Cc(k)^4 iv.kcneg; iv.kbpos*Cc(k)^3 -iv.kbneg 0 0; iv.kaneg 0 -iv.kapos*Cc(k)^4 0; iv.kcpos 0 0 -iv.kcneg];
            v = [vec(k); vec(k+nx); vec(k+2*nx); vec(k+3*nx)];
            v = linsolve((eye(4)-A*iv.dt),v);
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
