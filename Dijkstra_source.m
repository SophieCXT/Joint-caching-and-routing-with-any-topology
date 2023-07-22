function [spcost sp] = Dijkstra_source(G, s)
% Find s->all destination shortest paths.
% This is an implementation of the dijkstras algorithm, wich finds the 
% minimal cost path (sp) and the path length (spcost) between two nodes. Its supoussed to solve the problem on 
% possitive weighted instances.
% NOTE: zero entries in matrix G mean edges
% do not exist (i.e., inf-cost); o.w. the entry (i,j) means link weight
% from node i to node j (assumed to be positive).

% the inputs of the algorithm are:
% G: the adjacency matrix (G(i,j)==1 if link (i,j) exists; G(i,j)==0 otherwise);
% s: source node index;
% d: destination node index;

% Output:
% spcost: n-dimensional array of s->d distance (hop count) for each d 
% sp: n-dimensional cell array of s->d shortest path for each d; each path
% is represented by a vector of node indices (including s and d)

% For information about this algorithm visit:
% http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

% Author: Ting He, March 28, 2014.


n=size(G,1);
S(1:n) = 0;     %s, vector, set of visited vectors
dist(1:n) = inf;   % it stores the shortest distance between the source node and any other node;
prev(1:n) = n+1;    % Previous node, informs about the best previous node known to reach each  network node 

dist(s) = 0;

while sum(S)~=n    
    Q = find(S==0); % unvisited nodes
    [min_dist u]=min(dist(Q)); u = Q(u); 
    S(u)=1;
    if min_dist == inf % s is not connected to u
        break;
    end
    for i=1:n
        if G(u,i) ~= 0 && ~S(i)
            alt = dist(u) + G(u,i);
            if alt < dist(i)
                dist(i)= alt;
                prev(i)=u;
            end
        end
    end
end

spcost = zeros(1,n);
sp = cell(1,n);
for d=1:n % for each destination:    
    if dist(d) == inf
        sp{d} = []; % no s->d path available
    else
        sp{d} = [d];
        while sp{d}(1) ~= s
            if prev(sp{d}(1))<=n
                sp{d}=[prev(sp{d}(1)) sp{d}];
            else
                disp('Error: Dijkstra fails to record path');
            end
        end
    end
    spcost(d) = dist(d);
end