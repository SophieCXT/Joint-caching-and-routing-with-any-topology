% Special case 1 of online settings: unlimited link capacities (routing = RNR, integral content placement)
% Data is from youtube viewing history
% topology = 'Abovenet';
% topology = 'Level3';
topology = 'VerioUS';

addpath("matlab_bgl");
% >> cd C:\gurobi912\win64\matlab
% >> gurobi_setup

switch topology
    case 'Abovenet'
        C = 10; % total #items, following [Ioannidis18JSAC]
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
%         skewness = 1.2; % skewness is demands, with average rate of 1 per client [Ioannidis18JSAC]
%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes
        C_client = 1:4; 
        deg_client = 3; % maximum degree of clients
        n_client = 5;
        k_paths = 10; % [Ioannidis18JSAC]
        K_paths = [1 10 20 30];
    case 'Level3'
        C = 15; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]     
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;        
        deg_client = 5; %17; % maximum degree of clients
        n_client = 8;
        k_paths = 10; % [Ioannidis18JSAC]
        K_paths = [1 10 20 30];
    case 'VerioUS'
        C = 50; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)        
        skewness = 1.2; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]                
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;                 
        deg_client = 1; %17; % maximum degree of clients
        n_client = 12; % #clients (with degree up to deg_client); note: this includes caches (but excludes the server)
        k_paths = 10; % [Ioannidis18JSAC]
        K_paths = [1 10 20 30];
end

load(['data/' topology '.mat']); % A: sparse adjacency matrix with binary entries
A = full(A);
V = length(A); 
window_size = 2;
weeks = 18;
fontsize = 12; % tian, at least ten monte carlo runs.

% index,views,size
M = readmatrix('data/size_total_views_desc.txt');
size = M(1:C,3);
min_size = min(size);
size_init = round(size ./ min_size);

%% vary edge cache capacity:
% k_paths = 10;
% C_client = [1 2 3 4 5]; % cache capacity at clients/edge nodes

C_RNR = zeros(length(C_client),3); % C_RNR(i,1): routing cost in setting i under proposed solution; C_RNR(i,2): Ioannidis' solution; C_RNR(i,3): en-route caching along shortest paths to the server
for n_week = 4:weeks
    disp(' ')
    disp(['n_week = ' num2str(n_week) ':'])
    disp(' ')
[ G,server,client,request ] = ewma_generate_instance_unlimited_link_capacities( n_week,window_size,A,cmin,cmax,cmin_s,cmax_s,n_client,C );
% generate candidate paths for [Ioannidis18JSAC]:
netCostMatrix = G;
netCostMatrix(A==0) = inf; 
candidatePaths = cell(1,V);
% method 2: find a k shortest paths for each requester for a fixed k
% k_paths = 10; % [Ioannidis18JSAC]
for i = 1:length(client)
    [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), k_paths); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
    candidatePaths{client(i)} = shortestPaths;
    disp(['client ' num2str(i) ': found ' num2str(k_paths) ' paths with stretch <= ' num2str(max(totalCosts)/min(totalCosts))])
    if totalCosts(1) > min(totalCosts)+10^(-6)
        error(['client ' num2str(i) ': the first k-shortest path is not the shortest path']);
    end
end% candidatePaths{i}{j}: j-th shortest path from node i (a client) to server
candidatePaths1 = cell(1,V);
for i=1:length(client)
    candidatePaths1{client(i)} = candidatePaths{client(i)}(1);
end
for i_cv = 1:length(C_client)
    c_client = C_client(i_cv);
    disp(' ')
    disp(['c_client = ' num2str(c_client)])
    [ c_RNR, c_RNR_Ioannidis, c_RNR_SP ] = simulate_unlimited_link_capacities( A,C,c_v,c_client,k_paths,  G,server,client,request, candidatePaths, candidatePaths1);
    C_RNR(i_cv,1) = C_RNR(i_cv,1)+c_RNR;
    C_RNR(i_cv,2) = C_RNR(i_cv,2)+c_RNR_Ioannidis;
    C_RNR(i_cv,3) = C_RNR(i_cv,3)+c_RNR_SP;
end
end
C_RNR = C_RNR./(weeks-window_size+1); 

save(['data/online_settings/ewma_no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology '.mat'], 'C_client','C_RNR'); 

% load(['data/no_link_capacity_cachecapacity_' topology '.mat']);
figure;
bar(C_client, C_RNR);
legend('Alg. 1','k shortest paths','shortest path');
xlabel('cache capacity');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
saveas(gcf,['data/online_settings/ewma_no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% vary k_paths:
% c_client = 2; % following [Ioannidis18JSAC]
% K_paths = [1 10 20 30];

C_RNR = zeros(length(K_paths),2); % C_RNR(i,1): routing cost in setting i under proposed solution; C_RNR(i,2): Ioannidis' solution; 
for n_week = window_size:weeks%:weeks 
[ G,server,client,request ] = ewma_generate_instance_unlimited_link_capacities( n_week,window_size,A,cmin,cmax,cmin_s,cmax_s,n_client,C );
V = length(A);
c_cache = ones(V,1)*c_v;
c_cache(server) = C; 
c_cache(client) = c_client;
[ c_RNR, var_x, var_r ] = integral_caching_RNR( G, c_cache, request );
C_RNR(:,1) = C_RNR(:,1) + c_RNR;
% pre-compute k-shortest paths for speed:
netCostMatrix = G;
netCostMatrix(A==0) = inf;
candidatePaths = cell(1,V);
for i = 1:length(client)
    [shortestPaths, totalCosts] = kShortestPath(netCostMatrix, server, client(i), max(K_paths)); % shortestPaths: 1*k_paths cell array, shortestPaths{i}: node sequence for i-th shortest path; totalCosts: costs of the paths as 1*k_paths array
    candidatePaths{client(i)} = shortestPaths;
    disp(['client ' num2str(i) ': found ' num2str(max(K_paths)) ' paths with stretch <= ' num2str(max(totalCosts)/min(totalCosts))])
end
for i_cv = 1:length(K_paths)
    k_paths = K_paths(i_cv);
    disp(' ')
    disp(['k_paths = ' num2str(k_paths)])
    candidatePaths1 = cell(1,V);
    for i = 1:length(client)
        candidatePaths1{client(i)} = candidatePaths{client(i)}(1:min(k_paths,length(candidatePaths{client(i)})));
    end
    [c_RNR_Ioannidis ] = integral_caching_RNR_Ioannidis( G, c_cache, request, candidatePaths1, client );
    C_RNR(i_cv,2) = C_RNR(i_cv,2)+c_RNR_Ioannidis;
end
end
C_RNR = C_RNR./(weeks-window_size+1); %runs; 

save(['data/online_settings/ewma_no_link_capacity_candidatepaths_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology '.mat'], 'K_paths','C_RNR','c_client'); 

figure;
bar(K_paths, C_RNR);
legend('Alg. 1','k shortest paths')
xlabel('#candidate paths');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
saveas(gcf,['data/online_settings/ewma_no_link_capacity_candidatepaths_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)