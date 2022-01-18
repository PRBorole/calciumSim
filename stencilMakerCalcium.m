function stencil = stencilMakerCalcium(filename)

[A,id,pid,coord,r,subset]=readSWC(filename);
[M,nLst,~,~,nx,~,meanEdge,~,~,~]=getGraphStructure(filename,false,false);
meanEdge = meanEdge*10^-5;

% s,t are for making the edges
s = id(2:end); t = pid(2:end);

% this creates a graph
G = graph(s,t);

% this is the adjacency matrix
M = adjacency(G);
M = full(M);
stencil = M + eye(nx);

for i = 1:nx 
    ptCoord = coord(i,:);
    nNeigh = size(nLst{i},1);
    dx = zeros(1,nNeigh);
    for j = 1:nNeigh
        idx = nLst{i}(j);
        %ptnCoord = coord(idx,:);
        %dx(j) = sqrt(sum((ptCoord - ptnCoord).^2));
        dx(j) = meanEdge;
        stencil(i,idx) = 1/dx(j);
    end
    
    
    stencil(i,:) = 2*stencil(i,:)/sum(dx);
    stencil(i,i) = -stencil(i,i)*sum(1./dx);
end
end

