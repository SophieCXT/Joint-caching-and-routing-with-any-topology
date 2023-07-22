function [ cost, load, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, Vs, request, K )
%Integral routing and source selection under binary cache capacities
%(caller of Algorithm 2):
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_link: V*V matrix containing link capacities; c_link(i,j) is the capacity of
% link (i,j). 
% Vs: array of nodes storing the (entire) catalog (rest do not store
% anything).
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% K: design parameter (assumed to be a positive integer). 
% Output:
% cost: total routing cost 
% load: V*V matrix; load(u,v) is the load factor (flow/capacity) for link (u,v) if the link exists (can be greater than 1
% due to bicriteria approximation)
% paths: 1*n cell array; paths{i}: the path for request i as a node
% sequence.
% request_type: n*2 matrix; request_type(i,:) describes the item and the
% requester for request of type i (routed on path paths{i}). 
% cost_splittable: minimum cost of splittable flow (that does not exceed
% any link capacity). 
% load_splittable: maximum load factor under splittable flow (must be <=
% 1). 

% generate auxiliary topology:
[C,V] = size(request); % #items and #nodes
G = [G zeros(V,1); zeros(1,V+1)];
G(V+1,Vs) = 1; % set virtual link costs = 1 (as G(.,.)=0 is interpreted as no link); the final cost needs to be reduced by sum(sum(request))*1 to ignore this cost of virtual links.
c_link = [c_link zeros(V,1); zeros(1,V+1)];
c_link(V+1,Vs) = sum(sum(request)); % unlimited capacity for virtual links
s = V+1;
n = nnz(request); 
request_type = zeros(n,2); % i-th request: item request_index(i,1), requested by node request_index(i,2)
lambda = zeros(n,1);
d = zeros(n,1);
j = 0;
for i=1:C
    for v=1:V
        if request(i,v) > 0
            j = j + 1;
            request_type(j,:) = [i v];
            lambda(j) = request(i,v);
            d(j) = v;
        end
    end
end
% calling Algorithm 2 to solve MSUFP:
[ cost, load_vector, paths, links, cost_splittable, load_vector_splittable ] = bicriteria_MSUFP_gurobi( G, c_link, s, d, lambda, K );
% convert solution back to original topology:
cost = cost - sum(lambda); % to ignore the cost incurred on virtual links
cost_splittable = cost_splittable - sum(lambda); 
load_splittable = 0; 
load = zeros(V);
for e=1:length(links(:,1))
    if max(links(e,:)) <= V % if it is a real link
        load(links(e,1),links(e,2)) = load_vector(e); % for our MSUFP algorithm, record the load for every link
        load_splittable = max(load_splittable, load_vector_splittable(e)); % for splittable flow, only record the maximum load of the busiest link
    end
end
for p=1:length(paths)
    paths{p} = paths{p}(2:end); % remove the virtual source from each path
end


end

