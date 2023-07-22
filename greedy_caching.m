function [ cost, load_max, placement, flows, paths, request_type  ] = greedy_caching( G, c_cache, c_link, request, integral, heuristic, placement )
%Greedy content placement based on optimal/suboptimal routing (and source
%selection):
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_cache: V*1 vector of cache capacities; c_cache(v): capacity of cache at
% node v.
% c_link: V*V matrix containing link capacities; c_link(i,j) is the capacity of
% link (i,j). 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% integral: 1 if integral source selection & routing; 0 if fractional
% heuristic: in the case of integral source selection & routing, 1
% indicating that heuristic is used to solve the MMUFP problem ('LP relaxation with randomized rounding'); if 0, solve ILP
% placement: initial placement to guarantee feasibility of the routing
% problem. 
% Output:

[C,V] = size(request); % #items and #nodes
c_cache = reshape(c_cache,V,1);
[ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement, request, integral, heuristic );% initial cost %cost = inf; load_max = inf; 
iter = 1;
while any(sum(placement,2)<c_cache) % if still space at some cache
    placement_min = placement; 
    V_notfull = find(sum(placement,2)<c_cache);
    for i_v = 1:length(V_notfull)
        v = V_notfull(i_v);
        C_notcached = find(placement(v,:)==0);
        for i_i = 1:length(C_notcached)
            i = C_notcached(i_i);
            placement1 = placement; placement1(v,i) = 1; 
            [ cost1, load_max1, flows1, paths1, request_type1 ] = source_selection_routing_given_placement( G, c_link, placement1, request, integral, heuristic );
            if (load_max>1 && (load_max1<load_max || (load_max1<=load_max && cost1<cost) )) ... % if originally congested: reduce congestion first
                    || (load_max<=1 && load_max1<=1 && cost1<cost) % if not congested: reduce cost
                cost=cost1; load_max=load_max1; placement_min=placement1; flows=flows1; paths=paths1; request_type=request_type1;
            end 
        end
    end
    disp(['greedy_caching: finish iteration ' num2str(iter) ', cost = ' num2str(cost) ', congestion = ' num2str(load_max)])
    iter = iter + 1; 
    if ~any(placement_min>placement) % no further placement will help
        disp(['greedy: stopping early as no further placement of content helps'])
        break;
    end
    placement = placement_min; 
end        

end

