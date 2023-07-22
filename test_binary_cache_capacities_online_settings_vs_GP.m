% Special case 2: binary cache capacities (reduced to MSUFP)
topology = 'Abovenet';
% topology = 'TiscaliEurope';
% topology = 'Level3';
% topology = 'VerioUS';

% addpath("matlab_bgl");
% addpath("C:\gurobi912\win64\matlab\");
% >> cd C:\gurobi912\win64\matlab
% >> gurobi_setup

switch topology
    case 'Abovenet'
        C = 10; % 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        linkmin = 13715.796; % .5; % minimum link capacity (before augmentation)
        linkmax = 13715.796; %equals 15 Gbps % maximum link capacity (before augmentation)
        Linkmin = [3 4 5 6].*13715.796./3;%[2 3 4 5].*13715.796./3;%5=n_client, [.1 .2 .3 .4] [.2 .3 .4 .5]; %[.5 .75 1 1.25 ];
        Linkmin_Gbps = [15 20 25 30];%[10 15 20 25]; %in Gbps
        K_total = [2:10:162]; % optimal K = 12
        K = 160; %12; % optimized for ours        
        n_cache = 1; % #caches storing the entire catalog
        deg_client = 3; %17; % maximum degree of clients
        n_client = 5;
        cache_pre = 0; % not predefine the cache
       case 'Level3'
        C = 15; % 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]%
        n_client = 8;
        linkmin = .2./n_client; %.4; % minimum link capacity (before augmentation)
        linkmax = .2./n_client; %.4; % maximum link capacity (before augmentation)
        Linkmin = [.2 .3 .4 .5]./n_client; %[.5 .75 1 1.25 ];
        K_total = [2:2:20]; % optimal K = ?
        K = 18; % optimized for ours        
        n_cache = 1; % #caches storing the entire catalog
        deg_client = 5; %17; % maximum degree of clients 
        cache_pre = 0; % not predefine the cache
    case 'VerioUS'
        C = 50; % 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]
%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        linkmin = .2; % minimum link capacity (before augmentation)
        linkmax = .2; % maximum link capacity (before augmentation)
        Linkmin = [.1 .2 .3 .4]; %[.5 .75 1 1.25 ];
        K_total = [2:8:50]; % optimal K = ?
        K = 18; % optimized for ours        
        n_cache = 1; % #caches storing the entire catalog
        deg_client = 1; %17; % maximum degree of clients
        n_client = 12; % #clients (with degree up to deg_client); note: this includes caches (but excludes the server)
        cache_pre = 0; % not predefine the cache        
end

% RNR represents what was proposed in [Ioannidis18JSAC]. 
load(['data/' topology '.mat']); % A: sparse adjacency matrix with binary entries
A = full(A);
V = length(A); 
K0 = 2; % fixed for benchmark [Skutella02]
fontsize = 16; %used to be 12
window_size = 2;
weeks = 18;
predict_batch_size = 5; % one hour = one Monte Carlo run, predicting 5 hours at a time
predict_batch_num = 20; %20, num of iterations, GP_top30_videos_model_plus_5_pred_11_29_[1~20]
w1 = .85; 
w2 = .4;

total_num_request=34289494;

% %% compute the bound on link load:
% [ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities( A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_cache,deg_client, linkmin, linkmax);
% c_link_min = min(c_link(c_link>0)); % minimum link capacity
% bound = 2.^(1./K_total).*c_link_min + max(request(request>0)).*2.^(1./K_total)./2./(2.^(1./K_total)-1); % upper bound on absolute link load for the min-capacity link
% disp(['lambda_max = ' num2str(max(request(request>0))) ', c_min = ' num2str(min(c_link(c_link>0))) ])
% load_bound = bound./c_link_min % upper bound on the maximum link load factor
% %%
% K=5; 
% [ cost, load, paths, request_type, cost_splittable,~ ] = integral_routing_binary_cache( G, c_link, [server; cache]', request, K );
% disp(['proposed: cost = ' num2str(cost) ' (min cost of splittable flow: ' num2str(cost_splittable) '), max load = ' num2str(max(max(load))) ])
%% vary design parameter K:

Cost = zeros(length(K_total),2); % Cost(:,1): proposed, Cost(:,2): splittable flow
Load = zeros(length(K_total),3); % 1: proposed, 2: splittable, 3: upper bound
tic % tian
% for n_week = window_size:weeks
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
[ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities_online_settings( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre);
% compute bound:
c_link_min = min(c_link(c_link>0)); % minimum link capacity
bound = 2.^(1./K_total).*c_link_min + max(request(request>0)).*2.^(1./K_total)./2./(2.^(1./K_total)-1); % upper bound on absolute link load for the min-capacity link
load_bound = bound./c_link_min; % upper bound on the maximum link load factor
Load(:,3) = Load(:,3) + load_bound'; 
% % compute RNR:
% placement = zeros(V,C);
% placement([server; cache], :) = 1;
% [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request );
% disp(['RNR: cost = ' num2str(cost) ', max load = ' num2str(load_max) ])
% Cost(:,2) = Cost(:,2)+cost;
% Load(:,2) = Load(:,2)+load_max; 
for i_p = 1:length(K_total)
    K = K_total(i_p);
disp(' ')
disp(['lambda_max = ' num2str(max(request(request>0))) ', c_min = ' num2str(min(c_link(c_link>0))) ', K = ' num2str(K) ':'])
[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache( G, c_link, [server; cache]', request, K );
disp(['proposed: cost = ' num2str(cost) ' (min cost of splittable flow: ' num2str(cost_splittable) '), max load = ' num2str(max(max(load_unsplit))) ' (load of splittable flow: ' num2str(load_splittable) ')'])

Cost(i_p,1) = Cost(i_p,1) + cost;
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit)); 
end
Cost(:,2) = Cost(:,2) + cost_splittable; 
Load(:,2) = Load(:,2) + load_splittable; 
end
end
Cost = Cost./total_num_request;%(predict_batch_num*predict_batch_size*total_num_request*1000*3600);%runs;
Load = Load./(predict_batch_num*predict_batch_size);%runs;
time = toc./(predict_batch_num*predict_batch_size);%runs; % tian: execution time
is_save = 1;
if is_save
save(['data/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_' topology '.mat'],'K_total','Cost', 'Load', 'linkmin','linkmax','time'); 
end

figure;
% bar(K_total,Cost);
plot(K_total,Cost(:,1),'b-',...
    K_total,Cost(:,2),'r-','LineWidth',2);
legend('Alg. 2', ...'RNR', 
    'splittable flow')
xlabel('K');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
if is_save
saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_' topology '_cost'],'epsc')
end
figure;
% bar(K_total,Load);
plot(K_total,Load(:,1),'b-',...
    K_total,Load(:,2),'r-',...
    ...K_total,Load(:,3),'c--',...
    'LineWidth',2);
legend('Alg. 2',...'RNR',
    'splittable flow') %,...
    %'upper bound') % the upper bound is very loose for large K
xlabel('K');
ylabel('congestion');
set(gca, 'FontSize', fontsize);
if is_save
saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_' topology '_congestion'],'epsc')
end
load handel;
player = audioplayer(y, Fs);
play(player)
% end% for preindx

%% vary (minimum) link capacity:
% K = 12; % optimized for ours
% K0 = 10; % fixed for benchmark
% Linkmin = [.5 .75 1 1.25 1.5 ];
alpha = 1; % linkmax/linkmin = alpha --> all links have identical capacities before augmentation

Cost = zeros(length(Linkmin),4); % 1: K, 2: K0, 3: splittable flow, 4: RNR
Load = zeros(length(Linkmin),6); % 1: K, 2: K0, 3: splittable flow, 4: RNR, 5: K bound, 6: K0 bound
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
[ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities_online_settings( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre);
     placement = zeros(V,C);
    placement([server; cache], :) = 1;
for i_p = 1:length(Linkmin)
    linkmin = Linkmin(i_p);
    linkmax = alpha*linkmin; 
    c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link;
    c_link = max(c_link,demand_SP); % augment link capacities to ensure it is feasible to serve each request
    disp(' ')
    disp(['lambda_max = ' num2str(max(request(request>0))) ', c_min = ' num2str(min(c_link(c_link>0))) ', c_max = ' num2str(max(c_link(c_link>0))) ':'])
    % compute bound:
    c_link_min = min(c_link(c_link>0)); % minimum link capacity
    bound = 2.^(1./[K K0]).*c_link_min + max(request(request>0)).*2.^(1./[K K0])./2./(2.^(1./[K K0])-1); % upper bound on absolute link load for the min-capacity link
    Load(i_p,5:6) = Load(i_p,5:6) + bound./c_link_min; % upper bound on the maximum link load factor
    % compute RNR:
    [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request );
    disp(['RNR: cost = ' num2str(cost) ', max load = ' num2str(load_max) ])
    Cost(i_p,4) = Cost(i_p,4)+cost;
    Load(i_p,4) = Load(i_p,4)+load_max;
[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, [server; cache]', request, K );
disp(['K = ' num2str(K) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,1) = Cost(i_p,1) + cost;
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit));

[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, [server; cache]', request, K0 );
disp(['K = ' num2str(K0) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,2) = Cost(i_p,2) + cost;
Load(i_p,2) = Load(i_p,2) + max(max(load_unsplit));

Cost(i_p,3) = Cost(i_p,3) + cost_splittable; 
Load(i_p,3) = Load(i_p,3) + load_splittable;
end
end
end
Cost = Cost./total_num_request;%(predict_batch_num*predict_batch_size*total_num_request*1000*3600);%runs;
Load = Load./(predict_batch_num*predict_batch_size);%runs;

save(['data/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_raw_' topology '.mat'],'Linkmin','Cost', 'Load', 'K','K0'); 
% 
% figure;
% bar(Linkmin_Gbps,Cost(:,[3 1 2 4]));
% legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR')
% xlabel('physical link capacity (Gbps)');
% ylabel('routing cost/request (ms)');
% set(gca, 'FontSize', fontsize);
% % title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
% saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_cost'],'epsc')
% 
% 
% 
% figure;
% bar(Linkmin_Gbps,Load(:,[3 1 2 4 ]));
% %hold on;
% %plot([min(Linkmin) max(Linkmin)],[1 1],'r--','LineWidth',2);
% %hold off;
% legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
% xlabel('physical link capacity (Gbps)');
% ylabel('congestion');
% ylim([0.5 max(Load(:,2))+.1])
% set(gca, 'FontSize', fontsize);
% title(['Topology ', topology])
% saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_congestion'],'epsc')
% 
% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% gaussian_process vary design parameter K:

Cost = zeros(length(K_total),2); % Cost(:,1): proposed, Cost(:,2): splittable flow
Load = zeros(length(K_total),3); % 1: proposed, 2: splittable, 3: upper bound
tic % tian
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
[ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities_gaussian_process( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre);
% compute bound:
c_link_min = min(c_link(c_link>0)); % minimum link capacity
bound = 2.^(1./K_total).*c_link_min + max(request(request>0)).*2.^(1./K_total)./2./(2.^(1./K_total)-1); % upper bound on absolute link load for the min-capacity link
load_bound = bound./c_link_min; % upper bound on the maximum link load factor
Load(:,3) = Load(:,3) + load_bound'; 
% % compute RNR:
% placement = zeros(V,C);
% placement([server; cache], :) = 1;
% [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request );
% disp(['RNR: cost = ' num2str(cost) ', max load = ' num2str(load_max) ])
% Cost(:,2) = Cost(:,2)+cost;
% Load(:,2) = Load(:,2)+load_max; 
for i_p = 1:length(K_total)
    K = K_total(i_p);
disp(' ')
disp(['lambda_max = ' num2str(max(request(request>0))) ', c_min = ' num2str(min(c_link(c_link>0))) ', K = ' num2str(K) ':'])
[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, [server; cache]', request, K );
disp(['proposed: cost = ' num2str(cost) ' (min cost of splittable flow: ' num2str(cost_splittable) '), max load = ' num2str(max(max(load_unsplit))) ' (load of splittable flow: ' num2str(load_splittable) ')'])

Cost(i_p,1) = Cost(i_p,1) + cost;
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit)); 
end
Cost(:,2) = Cost(:,2) + cost_splittable; 
Load(:,2) = Load(:,2) + load_splittable; 
end
end
Cost = Cost./(predict_batch_num*predict_batch_size*total_num_request*1000*3600);%runs;
Load = Load./(predict_batch_num*predict_batch_size);%runs;
time = toc./(predict_batch_num*predict_batch_size);%runs; % tian: execution time
is_save = 1;
if is_save
save(['data/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_gaussian_process_' topology '.mat'],'K_total','Cost', 'Load', 'linkmin','linkmax','time'); 
end

figure;
% bar(K_total,Cost);
plot(K_total,Cost(:,1),'b-',...
    K_total,Cost(:,2),'r-','LineWidth',2);
legend('Alg. 2', ...'RNR', 
    'splittable flow')
xlabel('K');
ylabel('routing cost/request (ms)');
title(['gaussian process: Topology ', topology, ', predict batch size=' num2str(predict_batch_size), ', predict batch num=' num2str(predict_batch_num)])
set(gca, 'FontSize', fontsize);
if is_save
saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_cost_gaussian_process'],'epsc')
end
figure;
% bar(K_total,Load);
plot(K_total,Load(:,1),'b-',...
    K_total,Load(:,2),'r-',...
    ...K_total,Load(:,3),'c--',...
    'LineWidth',2);
legend('Alg. 2',...'RNR',
    'splittable flow') %,...
    %'upper bound') % the upper bound is very loose for large K
xlabel('K');
ylabel('congestion');
set(gca, 'FontSize', fontsize);
% title(['gaussian process: Topology ', topology, ', predict_batch_size=' num2str(predict_batch_size), ', predict_batch_num=' num2str(predict_batch_num)])
if is_save
saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin)  '_' topology '_congestion_gaussian_process'],'epsc')
end
load handel;
player = audioplayer(y, Fs);
play(player)
% end% for preindx

%% gaussian_process vary (minimum) link capacity:
% K = 12; % optimized for ours
% K0 = 10; % fixed for benchmark
% Linkmin = [.5 .75 1 1.25 1.5 ];
alpha = 1; % linkmax/linkmin = alpha --> all links have identical capacities before augmentation

Cost = zeros(length(Linkmin),4); % 1: K, 2: K0, 3: splittable flow, 4: RNR
Load = zeros(length(Linkmin),6); % 1: K, 2: K0, 3: splittable flow, 4: RNR, 5: K bound, 6: K0 bound

% n_timestamp=0;
% batch_num=0;
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
[ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities_gaussian_process( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre);
    placement = zeros(V,C);
    placement([server; cache], :) = 1;
for i_p = 1:length(Linkmin)
    linkmin = Linkmin(i_p);
    linkmax = alpha*linkmin; 
    c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link;
    c_link = max(c_link,demand_SP); % augment link capacities to ensure it is feasible to serve each request
    disp(' ')
    disp(['lambda_max = ' num2str(max(request(request>0))) ', c_min = ' num2str(min(c_link(c_link>0))) ', c_max = ' num2str(max(c_link(c_link>0))) ':'])
    % compute bound:
    c_link_min = min(c_link(c_link>0)); % minimum link capacity
    bound = 2.^(1./[K K0]).*c_link_min + max(request(request>0)).*2.^(1./[K K0])./2./(2.^(1./[K K0])-1); % upper bound on absolute link load for the min-capacity link
    Load(i_p,5:6) = Load(i_p,5:6) + bound./c_link_min; % upper bound on the maximum link load factor
    % compute RNR:
    [ cost, load_max ] = RNR_given_placement( G, c_link, placement, request );
    disp(['RNR: cost = ' num2str(cost) ', max load = ' num2str(load_max) ])
    Cost(i_p,4) = Cost(i_p,4)+cost;
    Load(i_p,4) = Load(i_p,4)+load_max;
[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, [server; cache]', request, K );
disp(['K = ' num2str(K) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,1) = Cost(i_p,1) + cost;
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit));

[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache_gurobi( G, c_link, [server; cache]', request, K0 );
disp(['K = ' num2str(K0) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,2) = Cost(i_p,2) + cost;
Load(i_p,2) = Load(i_p,2) + max(max(load_unsplit));

Cost(i_p,3) = Cost(i_p,3) + cost_splittable; 
Load(i_p,3) = Load(i_p,3) + load_splittable;
end
end
end
Cost_gaussian_process = Cost./total_num_request;%(predict_batch_num*predict_batch_size*total_num_request*1000*3600);%runs;
Load_gaussian_process = Load./(predict_batch_num*predict_batch_size);%runs;

save(['data/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process.mat'],'Linkmin','Cost_gaussian_process', 'Load_gaussian_process', 'K','K0'); 


% figure;
% bar(Linkmin_Gbps,Cost_gaussian_process(:,[3 1 2 4]));
% legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR')
% xlabel('physical link capacity (Gbps)');
% ylabel('routing cost/request (ms)');
% set(gca, 'FontSize', fontsize);
% % title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
% saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_' topology '_cost_gaussian_process'],'epsc')
% figure;
% bar(Linkmin_Gbps,Load_gaussian_process(:,[3 1 2 4 ]));
% %hold on;
% %plot([min(Linkmin) max(Linkmin)],[1 1],'r--','LineWidth',2);
% %hold off;
% legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
% xlabel('physical link capacity (Gbps)');
% ylabel('congestion');
% ylim([0.5 max(Load_gaussian_process(:,2))+.1])
% set(gca, 'FontSize', fontsize);
% title(['Topology ', topology])
% saveas(gcf,['plot/100_runs_trace_driven/case2_binary_cache_capacity_linkcapacity_K_' num2str(K) '_C_' num2str(C) '_' topology '_congestion_gaussian_process'],'epsc')
% 
% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% plot gaussian_process vs online settings
% running sequence: first part, then gaussian_process vary (minimum) link capacity,
% then first part, then 'vary (minimum) link capacity'
figure;
b=bar(Linkmin_Gbps,Load(:,[3 1 2 4 ]),w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
hold on
%gaussian_process
b=bar(Linkmin_Gbps,Load_gaussian_process(:,[3 1 2 4 ]),w2,'FaceAlpha',1);%bar(C_client, C_RNR,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
grid on
xlabel('physical link capacity (Gbps)');
ylabel('congestion');
ylim([0.5 max(Load(:,2))+.1])
set(gca, 'FontSize', fontsize);
saveas(gcf,['data/100_runs_trace_driven/solid_case2_congestion_gaussian_process_vs_online_binary'],'epsc')


figure;
b=bar(Linkmin_Gbps,Cost(:,[3 1 2 4 ]),w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
hold on
%gaussian_process
b=bar(Linkmin_Gbps,Cost_gaussian_process(:,[3 1 2 4 ]),w2,'FaceAlpha',1);%bar(C_client, C_RNR,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green, purple:[0.4940, 0.1840, 0.5560]
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
grid on
xlabel('physical link capacity (Gbps)');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
saveas(gcf,['data/100_runs_trace_driven/solid_case2_routing_cost_gaussian_process_vs_online_binary'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)
