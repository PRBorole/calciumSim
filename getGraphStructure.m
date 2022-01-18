function [M,nLst,bLst,brchLst,numNodes,numEdges, meanEdge,maxEdge,minEdge,medEdge]=getGraphStructure(filename, plt,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will construct the graph structure
% and output edge information, neigbor list
% boundary nodes and more
%-------------------------------------------------------------------------%
% Input: filename of SWC file, and true/false for plotting cell geometry
% and HINES sparsity pattern, verbose show printed output true/false
%
% Output: M is the adjacency matrix of the graph
%         nLst is the neighbor list for each node
%         bLst is the list of boundary nodes
%         brchLst is the list of branching nodes
%         numNodes, numEdges are number of nodes, number of edges
%         meanEdge, maxEdge, minEdge, medEdge are the
%         mean, maximum,minimum, and median edge lengths
%-------------------------------------------------------------------------%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the SWC file
% only need id, pid for graph structure
% and coord for ploting
[~,id,pid,coord,~,~]=readSWC(filename);

% s,t are for making the edges
s = id(2:end); t = pid(2:end);

% this creates a graph
G = graph(s,t);
% this is the adjacency matrix
M = adjacency(G);

% if plotting is turned on then plot
if plt
    % this makes a plot of the graph
    figure
    plot(G,'XData',coord(1:end,1),'YData',coord(1:end,2),'ZData',coord(1:end,3)); 
    
    % plot sparsity pattern of adjacency matrix
    % I add the identity because the adjacency matrix has zero main
    % diagonal, this matrix is the HINES matrix for the graph
    figure
    spy(M+eye(size(M)))
end

% next we get the neigborlists at every node
% use the 
nLst={};    % neighbor list for non boundary, non branch points
bLst={};    % boundary list
brchLst={}; % branch node list
    for i=1:height(G.Nodes)
        % append neighbor arrays to list
        nLst{end+1}=neighbors(G,i);
        if verbose
            fprintf('Node %i has nLst size = %i neighbor(s)\n',i,length(neighbors(G,i)))
        end
        % if the length of neighbor array is one then that node is a 
        % boundary node
        if length(neighbors(G,i))==1
            bLst{end+1}=i;
        end
        % if a node has more than two neighbors
        % then it is a branch point
        if length(neighbors(G,i))>2
            brchLst{end+1}=i;
        end
    end
    
% next we need cell structure
% number of nodes, number of edges, average edge length, max edge length..
numEdges = height(G.Edges); numNodes = height(G.Nodes);

% initialize array of edge length
edge_lengths=zeros(1,numEdges);
    for i=2:numNodes
        % I use the id and pid to find the length of each edge
        diff = (coord(id(i),:)-coord(pid(i),:)).^2;
        edge_lengths(i-1) = sum(diff)^0.5;
        if verbose
            fprintf('The edge (%i,%i) has length = %f  microns!\n',id(i),pid(i),edge_lengths(i-1))
    
        end
    end
% edge data
meanEdge = mean(edge_lengths); maxEdge = max(edge_lengths);
minEdge = min(edge_lengths); medEdge = median(edge_lengths);
if verbose
    fprintf ('\n\n The aver edge length = %f microns \n',meanEdge)
    fprintf ('The max edge length = %f microns \n',maxEdge)
    fprintf ('The min edge length = %f microns \n',minEdge)
    fprintf ('The median edge length = %f microns \n',medEdge)
end
end

