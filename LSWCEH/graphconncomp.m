function [s,c]=graphconncomp(G,NVPArgs)
%GRAPHCONNCOMP finds the connected components in graph.
% 
% [S,C] = GRAPHCONNCOMP(G) finds the strongly connected components using
% Tarjan's algorithm. A strongly connected component is a maximal group of
% nodes that are mutually reachable without violating the edge directions.
% G is a square adjacency matrix that represents the directed graph; all
% nonzero entries indicate the presence of an edge. S is the number of
% components and C is a vector indicating to which component each node
% belongs.
% 
% Tarjan's algorithm time complexity is O(n+e), where n and e are number of
% nodes and edges respectively.  
% 
% GRAPHCONNCOMP(...,'DIRECTED',false) assumes G is an undirected graph, the
% upper triangle of the sparse matrix G is ignored. A DFS-based algorithm
% computes the connected components. Time complexity is O(n+e).
% 
% GRAPHCONNCOMP(...,'WEAK',true) finds the weakly connected components. A
% weakly connected component is a maximal group of nodes which are mutually
% reachable by violating the edge directions. Default is false. The state
% of this parameter has no effect on undirected graphs. Time complexity is
% O(n+e).
% 
% Remarks: 
%   - Note that by definition a single node can also be a strongly
%   connected component. 
%   - A DAG must not have any strongly connected components larger than
%   one. 
% 
% Example:
%    % Create a directed graph with 10 nodes and 17 edges
%    g = sparse([1 1 1 2 2 3 3 4 5 6 7 7 8 9 9  9 9], ...
%               [2 6 8 3 1 4 2 5 4 7 6 4 9 8 10 5 3],true,10,10)
%    h = view(biograph(g));
%    % Find the strongly connected components
%    [S,C] = graphconncomp(g)
%    % Mark the nodes for each component with different color
%    colors = jet(S);
%    for i = 1:numel(h.nodes)
%        h.Nodes(i).Color = colors(C(i),:);
%    end
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHISDAG, GRAPHISOMORPHISM,
% GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTOPOORDER, GRAPHTRAVERSE.
%
% References: 
%  [1]	R. E. Tarjan "Depth first search and linear graph algorithms" SIAM
%       Journal on Computing, 1(2):146-160, 1972. 
%  [2]  R. Sedgewick "Algorithms in C++, Part 5 Graph Algorithms"
%       Addison-Wesley, 2002. 

%   Copyright 2006-2021 The MathWorks, Inc.

arguments
    G {bioinfo.internal.ValidationHelper.mustBeAdjacencyGraph}
    NVPArgs.Directed(1,1) logical = true;
    NVPArgs.Weak(1,1) logical = false;
end

G_graph = bioinfo.internal.biograph2matlab(G, 'Directed', NVPArgs.Directed);

compArgs = {};
if NVPArgs.Weak
    compArgs = {'Type', 'weak'};
end

c = conncomp(G_graph, compArgs{:});
if isempty(c)
    s = 0;
else
    s = max(c);
end