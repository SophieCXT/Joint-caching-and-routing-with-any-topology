function [ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths )
% simulate algorithms applicable to the general case with limited cache capacities and limited link capacities:
% IC+IR:
disp('IC+IR:')
integral = 1; % integral routing (i.e., single-path routing)?
heuristic = 1; % compute routing by heuristic (only if integral routing)?
[C,V] = size(request);
% 1.1 proposed alternating optimization heuristic:
placement = zeros(V,C); placement(server,:) = 1; % initially only place contents at the remote server
cost_final = inf; load_final = inf; 
tic
for i=1:max_iter
    disp(['iteration ' num2str(i) ':'])
    [ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement, request, integral, heuristic );
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing routing'])
    [ cost, load_max, placement ] = integral_caching_given_paths( G, c_cache, c_link, request, flows, paths, request_type ); 
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing caching'])
    disp(' ')
    if cost >= cost_final && load_max >= load_final
        disp(['no improvement after ' num2str(i) ' iterations; stop'])
        break;
    else% improve at least one of (cost, load)
        if i==2 ... % first time optimizing content placement
                || (load_final > 1.01 && load_max < load_final - .01 && cost < cost_final*1.1)... if was congested and can reduce congestion (without increasing cost too much)
                || ( (load_max <= 1.01 || load_max <= load_final) && cost < cost_final  )... % or not congested and can reduce cost
                || ( load_max < load_final + .1 && cost < cost_final*.75 ) % or not much more congested and can significantly reduce cost
            cost_final = cost; load_final = load_max;
        end
    end
end
iter_alternating = min(i,max_iter); 
disp(['alternating optimization takes ' num2str(toc) ' sec'])
disp(' ')
cost_alternating = cost_final; load_alternating = load_final; 

% % greedy heuristic (greedy content placement + routing by 'source_selection_routing_given_placement'):
% placement = zeros(V,C);
% placement(server,:) = 1;
% tic
% [ cost, load_max, placement, flows, paths, request_type  ] = greedy_caching( G, c_cache, c_link, request, integral, heuristic, placement );
% disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after greedy caching'])
% disp(['greedy optimization takes ' num2str(toc) ' sec']) % about 80 times
% slower than alternating optimization, hence skipped

% 1.2 [Ioannidis'18TON] (shortest path -> content placement):
[SPcost, SP] = Dijkstra_source(G, server); 
n = nnz(request); % #commodities (types of requests with nonzero rate)
request_type = zeros(n,2); % request_type(i,:) describes type-i request
flows = zeros(1,n);
paths = cell(1,n);
j = 0;
for i=1:C
    for v=1:V
        if request(i,v) > 0
            j = j + 1;
            request_type(j,:) = [i v];
            flows(j) = request(i,v);
            paths{j} = SP{v};
        end
    end
end
[ cost, load_max, placement ] = integral_caching_given_paths( G, c_cache, c_link, request, flows, paths, request_type );
cost_SP = cost; load_SP = load_max;
placement_SP = placement; % record the content placement for enroute-caching along shortest path

% 1.3 [Ioannidis'18TON]+MMUFP (shortest path -> content placement -> MMUFP
% heuristic):
[ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement_SP, request, integral, heuristic );
cost_SP1 = cost; load_SP1 = load_max; 

% 1.4 [Ioannidis'18TON]+RNR (shortest path -> content placement -> RNR)
[ cost, load_max ] = RNR_given_placement( G, c_link, placement_SP, request );
cost_SP2 = cost; load_SP2 = load_max; 
% disp(['RNR: cost = ' num2str(cost) ', max load = ' num2str(load_max) ])

% 1.5 [Ioannidis'18JSAC] (k-shortest paths -> content placement -> RNR)
netCostMatrix = G;
netCostMatrix(G==0) = inf; 
candidatePaths = cell(1,V);
for i = 1:length(client)
    [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), k_paths); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
    candidatePaths{client(i)} = shortestPaths;    
end% candidatePaths{i}{j}: j-th shortest path from node i (a client) to server
[c_RNR_Ioannidis, placement, var_r ] = integral_caching_RNR_Ioannidis( G, c_cache, request, candidatePaths, client );
% [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request );
load = zeros(V); 
for i=1:C
    for s=1:length(client)
        p = find(var_r(:,s,i)>0); % candidatePaths{client(s)}{p} is the selected path to serve request(i,client(s))
        for k=1:length(candidatePaths{client(s)}{p})-1
            load(candidatePaths{client(s)}{p}(k),candidatePaths{client(s)}{p}(k+1)) = load(candidatePaths{client(s)}{p}(k),candidatePaths{client(s)}{p}(k+1)) + request(i,client(s))*prod(1-placement(candidatePaths{client(s)}{p}(k+1:end),i)) / c_link(candidatePaths{client(s)}{p}(k),candidatePaths{client(s)}{p}(k+1)); 
        end
    end
end
cost_kSP_RNR = c_RNR_Ioannidis; load_kSP_RNR = max(max(load)); 

% IC+FR:
disp('IC+FR:')
integral = 0;
% 2.1 proposed alternating optimization heuristic:
placement = zeros(V,C); placement(server,:) = 1; % initially only place contents at the remote server
cost_final = inf; load_final = inf; 
tic
for i=1:max_iter
    disp(['iteration ' num2str(i) ':'])
    [ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement, request, integral, heuristic );
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing routing'])
    [ cost, load_max, placement ] = integral_caching_given_paths( G, c_cache, c_link, request, flows, paths, request_type ); 
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing caching'])
    disp(' ')
    if cost >= cost_final
        disp(['no improvement after ' num2str(i) ' iterations; stop'])
        break;
    else
        cost_final = cost; load_final = load_max; 
    end
end
iter_alternating_FR = min(i,max_iter); 
disp(['alternating optimization takes ' num2str(toc) ' sec'])
disp(' ')
cost_alternating_FR = cost_final; load_alternating_FR = load_final; 
% 2.2 [Ioannidis'18TON]+MMUFP (shortest path -> content placement -> MMSFP:
[ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement_SP, request, integral, heuristic );
cost_SP1_FR = cost; load_SP1_FR = load_max;



end

