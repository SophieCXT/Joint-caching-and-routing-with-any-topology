function [ c_RNR, c_RNR_Ioannidis, c_RNR_SP ] = simulate_unlimited_link_capacities( A,C,c_v,c_client,k_paths, G,server,client,request, candidatePaths, candidatePaths1)
V = length(A);

c_cache = ones(V,1)*c_v;
c_cache(server) = C; 
c_cache(client) = c_client;
% %% generate candidate paths for [Ioannidis18JSAC]:
% netCostMatrix = G;
% netCostMatrix(A==0) = inf; 
% candidatePaths = cell(1,V);
% % % method 1: find all paths within given stretch (too slow)
% % stretch = 4; 
% % for i = 1:length(client)
% %     disp(['client ' num2str(i) ':'])
% %     k_paths = 1;
% %     while true
% %         [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), k_paths); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
% %         if max(totalCosts)/min(totalCosts) <= stretch
% %             disp(['found ' num2str(k_paths) ' paths with stretch <= ' num2str(max(totalCosts)/min(totalCosts))])
% %             candidatePaths{i} = shortestPaths;
% %             k_paths = k_paths+1;
% %         else
% %             break;
% %         end
% %     end
% % end
% % method 2: find a k shortest paths for each requester for a fixed k
% % k_paths = 10; % [Ioannidis18JSAC]
% for i = 1:length(client)
%     [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), k_paths); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
%     candidatePaths{client(i)} = shortestPaths;
%     disp(['client ' num2str(i) ': found ' num2str(k_paths) ' paths with stretch <= ' num2str(max(totalCosts)/min(totalCosts))])
%     if totalCosts(1) > min(totalCosts)+10^(-6)
%         error(['client ' num2str(i) ': the first k-shortest path is not the shortest path']);
%     end
% end% candidatePaths{i}{j}: j-th shortest path from node i (a client) to server
% candidatePaths1 = cell(1,V);
% for i=1:length(client)
%     candidatePaths1{client(i)} = candidatePaths{client(i)}(1);
% end
%%
tic
[ c_RNR, var_x, var_r ] = integral_caching_RNR( G, c_cache, request );
disp(['integral_caching_RNR finishes in ' num2str(toc) ' seconds'])
tic
[c_RNR_Ioannidis ] = integral_caching_RNR_Ioannidis( G, c_cache, request, candidatePaths, client );
disp(['integral_caching_RNR_Ioannidis(' num2str(k_paths) ' candidate paths) finishes in ' num2str(toc) ' seconds'])
tic
[c_RNR_SP ] = integral_caching_RNR_Ioannidis( G, c_cache, request, candidatePaths1, client ); % a request can only be served from caches along the shortest path to the server
disp(['integral_caching_RNR_Ioannidis(1 candidate path) finishes in ' num2str(toc) ' seconds'])

disp(['c_RNR = ' num2str(c_RNR) ', c_RNR_Ioannidis = ' num2str(c_RNR_Ioannidis) ', c_RNR_SP = ' num2str(c_RNR_SP)]);


end

