function [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general( A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths)
% generate a random instance of the problem (link costs, link capacities,
% request rates):
V = length(A);
% generate link costs:
G = zeros(size(A));
G(A>0) = cmin + (cmax-cmin).*rand(nnz(A),1);
% select client/server:
degree = sum(A,2); 
[~,I] = sort(degree); 
server = I(1); % server is the first node with the minimum degree %find(degree<=deg_server,1,'first'); 
% client = setdiff(find(degree<=deg_client),server);
client = I(2:n_client+1);
G(server,G(server,:)>0) = cmin_s + (cmax_s-cmin_s).*rand(1,sum(G(server,:)>0)); % server is remote, connected to test of network by links of delays in [100, 200]ms
% generate demands:
request_method = 3;
switch request_method
    case 1
        % method 1: Zipf over all C*#client request types
        request1 = [1:C*length(client)].^(-skewness);
        request1 = length(client).* request1 ./ sum(request1);
        request1 = request1(randperm(C*length(client))); % randomly permuted lambda_{i,s} with total rate = #clients
        request = zeros(C,V);
        request(:,client) = reshape(request1,C,length(client));
    case 2
        % method 2: Zipf over the C contents, equal for all clients
        request1 = [1:C].^(-skewness);
        request1 = length(client).* request1 ./ sum(request1); % request1(i): request rate for content i over all the clients
        request = zeros(C,V);
        for i=1:C
            request(i,client) = (request1(i)/length(client))*ones(1,length(client));
        end
    case 3
        % method 3: Zipf over the C contents at each client, with randomly
        % permuted ranking for each client
        request1 = [1:C].^(-skewness);
        request1 = request1' ./ sum(request1); % rates for requesting contents 1,...,C at each client (sum to one)
        request = zeros(C,V);
        for s=1:length(client)
            request(:,client(s)) = request1(randperm(C));
        end  
end
% assign cache capacities:
c_cache = ones(V,1)*c_v;
c_cache(server) = C; 
c_cache(client) = c_client;
% generate link capacities:
c_link = zeros(size(A));
linkmin = max(linkmin, max(request1)); % make sure that minimum link capacity >= maximum demand as typically assumed in the literature of MSUFP
rand_c_link = rand(nnz(A),1); % store the random numbers drawn for computing link capacities (before augmentation) to increase linkmax proportionally
c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link; %rand(nnz(A),1);

% augment link capacities to ensure it is feasible to serve each request
% 1. from the server along the shortest path:
[~,SP] = Dijkstra_source(G, server);
% 2. a randomly selected cycle-free path:
netCostMatrix = G;
netCostMatrix(A==0) = inf;
SP = cell(1,V);
for i=1:length(client)
    [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), k_paths); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
%     if max(totalCosts) ~= totalCosts(end)
%         error(['The k-th shortest path is not the longest among k-shortest paths']);
%     end
    SP{client(i)} = shortestPaths{end}; % augment the 'k_paths'-th shortest path from each client to the server
end

demand_SP = zeros(size(A)); % demand_SP(u,v): total load on link (u,v) if serving every request from the server along the shortest path
for s=1:length(client)
    demand = sum(request(:,client(s)));
    for k=1:length(SP{client(s)})-1
        demand_SP(SP{client(s)}(k),SP{client(s)}(k+1)) = demand_SP(SP{client(s)}(k),SP{client(s)}(k+1)) + demand;
    end
end
c_link = max(c_link,demand_SP);

if false
%% plot the topology to find out where to put the server and the clients:
smin = 4; smax = 10;
g = graph(A);
figure;
h = plot(g, 'Layout', 'force');
highlight(h,server,'NodeColor','g','MarkerSize',smin+(smax-smin)*1); % green for server
% highlight(h,cache,'NodeColor','c','MarkerSize',smin+(smax-smin)*.5); % cyan for caches
highlight(h,client,'NodeColor','r','MarkerSize',smin+(smax-smin)*.5); % red for clients
switches = setdiff(1:V,[server; client]); 
highlight(h,switches,'NodeColor','b', 'MarkerSize',smin); % blue for switches      
c_link_max = max(c_link(c_link>0));
for s=1:V
    for t=1:V
        if A(s,t)>0
            highlight(h,s,t,'LineWidth', .5+ 2*c_link(s,t)/c_link_max);% width (on top of 0.5) is proportional to link capacity
        end
        if demand_SP(s,t)>0
            highlight(h,s,t,'EdgeColor','r'); % red links are on the shortest paths from clients to server
        end
    end
end    

end


end

