function [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request )
%Source selection and routing under a given integral content placement by "route to the nearest replica" (ignoring link capacities):
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_link: V*V matrix containing link capacities; c_link(i,j) is the capacity of
% link (i,j). 
% placement: V*C matrix; placement(v,i) = x_{vi} (whether node v stores
% item i). 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% Output:
% cost: total routing cost of the solution
% load_max: load factor (load/capacity) for the most congested link

% generate auxiliary topology:
[C,V] = size(request); % #items and #nodes
G = [G zeros(V,C); zeros(C,V+C)];
c_link = [c_link zeros(V,C); zeros(C,V+C)];
total_request = sum(sum(request));
for i=1:C
    G(V+i,placement(:,i)==1) = 1; % virtual link (v_i, v) for each v storing i, with cost 1 (to be removed later)
    c_link(V+i,placement(:,i)==1) = total_request; % unlimited capacity on virtual links
end% V+i is the virtual source for content i
% RNR:
cost = 0; 
link_load = zeros(size(G));
for i=1:C
    [SPcost,SP] = Dijkstra_source(G, V+i);
    for s=1:V
        if request(i,s) > 0
            cost = cost + request(i,s)*SPcost(s);
            for k=1:length(SP{s})-1
                link_load(SP{s}(k),SP{s}(k+1)) = link_load(SP{s}(k),SP{s}(k+1)) + request(i,s)/c_link(SP{s}(k),SP{s}(k+1));
            end
        end
    end
end
% convert solution back to original topology:
cost = cost - sum(request(request>0));  % ignore the cost on virtual links
link_load = link_load(1:V,1:V);
load_max = max(max(link_load)); 


end

