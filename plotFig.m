function plotFig(tSol,ySol,nx,dt,list)

if ismember("Wave",list)
    % wave
    figure;
    plot(linspace(0,64,nx),ySol(1500/dt:10:4000/dt,1:nx))
    title("Cc");
    % legend();
end

if ismember("Cc",list)
    figure;
    plot(tSol,ySol(:,1:nx))
    title("Cc");
    % legend();
end

if ismember("b",list)
    figure;
    plot(tSol,ySol(:,nx+1:2*nx))
    title("b");
    % legend();
end

if ismember("Ce",list)
    figure;
    plot(tSol,ySol(:,2*nx+1:3*nx))
    title("Ce");
    % legend();
end

if ismember("p",list)
    figure;
    plot(tSol,ySol(:,3*nx+1:4*nx))
    title("p");
    % legend();
end

if ismember("be",list)
    figure;
    plot(tSol,ySol(:,4*nx+1:5*nx))
    title("be");
    % legend();
end

end