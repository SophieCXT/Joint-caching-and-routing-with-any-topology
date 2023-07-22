% Special case 2: binary cache capacities (reduced to MSUFP)
topology = 'Abovenet';
% topology = 'TiscaliEurope';
% topology = 'Level3';
% topology = 'VerioUS';

%addpath("matlab_bgl");

switch topology
    case 'Abovenet'
        C = 10; % 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        linkmin = .2; % .5; % minimum link capacity (before augmentation)
        linkmax = .2; %.5; % maximum link capacity (before augmentation)
        Linkmin = [.1 .2 .3 .4]; %[0.2 0.3 0.4 0.5][.5 .75 1 1.25 ];
        Linkmin_Gbps = [7.5 15 22.5 30];
        K_total = [2:2:20]; % optimal K = 12
        K = 18; %12; % optimized for ours        
        n_cache = 1; % #caches storing the entire catalog
        deg_client = 3; %17; % maximum degree of clients
        n_client = 5;
        cache_pre = 0; % not predefine the cache
       case 'Level3'
        C = 15; % 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]
%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        linkmin = .2; %.4; % minimum link capacity (before augmentation)
        linkmax = .2; %.4; % maximum link capacity (before augmentation)
        Linkmin = [.2 .3 .4 .5]; %[.5 .75 1 1.25 ];
        K_total = [2:2:20]; % optimal K = ?
        K = 18; % optimized for ours        
        n_cache = 1; % #caches storing the entire catalog
        deg_client = 5; %17; % maximum degree of clients
        n_client = 8; 
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
        K_total = [2:2:20]; % optimal K = ?
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
runs = 100; % how many times for monte carlo run
fontsize = 16;
face_alpha=0.8;
w1 = .85; 
w2 = .4;

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
% linkmin = .5; % minimum link capacity (before augmentation)
% linkmax = .5; % maximum link capacity (before augmentation)
% K_total = [2:2:20]; % optimal K = 12
 
% Cache_pre = setdiff(find(sum(A,2)==1),[server 27 34 36 38 48]);
% for preindx = 1:length(Cache_pre)
%     cache_pre = Cache_pre(preindx);
%     disp(' ')
%     disp(['cache = node' num2str(cache_pre) ':'])


Cost = zeros(length(K_total),2); % Cost(:,1): proposed, Cost(:,2): splittable flow
Load = zeros(length(K_total),3); % 1: proposed, 2: splittable, 3: upper bound
tic % tian
for rindx = 1:runs
[ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities( A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre);
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

Cost(i_p,1) = Cost(i_p,1) + cost./sum(sum(request));
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit)); 
end
Cost(:,2) = Cost(:,2) + cost_splittable./sum(sum(request)); 
Load(:,2) = Load(:,2) + load_splittable; 
end
Cost = Cost./runs;
Load = Load./runs; 
time = toc./runs; % tian: execution time
is_save = 1;
if is_save
save(['data/100_runs/binary_cache_capacity_K_C_' num2str(C) '_' topology '.mat'],'K_total','Cost', 'Load', 'linkmin','linkmax','time'); 
end

figure;
% bar(K_total,Cost);
plot(K_total,Cost(:,1),'b-',...
    K_total,Cost(:,2),'r-','LineWidth',2);
legend('Alg. 2', ...'RNR', 
    'splittable flow')
xlabel('K');
ylabel('routing cost/request (ms)');
%title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
if is_save
saveas(gcf,['plot/100_runs/binary_cache_capacity_K_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_cost'],'epsc')
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
%title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
if is_save
saveas(gcf,['plot/100_runs/binary_cache_capacity_K_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_congestion'],'epsc')
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
for rindx = 1:runs    
    [ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities( A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, min(Linkmin), max(Linkmin), cache_pre);
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
    Cost(i_p,4) = Cost(i_p,4)+cost./sum(sum(request));
    Load(i_p,4) = Load(i_p,4)+load_max;
[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache( G, c_link, [server; cache]', request, K );
disp(['K = ' num2str(K) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,1) = Cost(i_p,1) + cost./sum(sum(request));
Load(i_p,1) = Load(i_p,1) + max(max(load_unsplit));

[ cost, load_unsplit, paths, request_type, cost_splittable, load_splittable ] = integral_routing_binary_cache( G, c_link, [server; cache]', request, K0 );
disp(['K = ' num2str(K0) ': cost = ' num2str(cost) ', max load = ' num2str(max(max(load_unsplit))) ])
Cost(i_p,2) = Cost(i_p,2) + cost./sum(sum(request));
Load(i_p,2) = Load(i_p,2) + max(max(load_unsplit));

Cost(i_p,3) = Cost(i_p,3) + cost_splittable./sum(sum(request)); 
Load(i_p,3) = Load(i_p,3) + load_splittable;
end
end
Cost = Cost./runs;
Load = Load./runs;

save(['data/100_runs/binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '.mat'],'Linkmin','Cost', 'Load', 'K','K0'); 


%load(['data/binary_cache_capacity_linkcapacity_' topology '.mat']);

figure;
bar(Linkmin_Gbps,Cost(:,[3 1 2 4]),w1,'FaceAlpha',face_alpha);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR')
xlabel('link capacity');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
saveas(gcf,['plot/100_runs/binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_cost'],'epsc')

figure;
bar(Linkmin_Gbps,Load(:,[3 1 2 4 ]),w1,'FaceAlpha',face_alpha);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
%hold on;
%plot([min(Linkmin) max(Linkmin)],[1 1],'r--','LineWidth',2);
%hold off;
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
xlabel('link capacity');
ylabel('congestion');
ylim([0.5 max(Load(:,2))+.1])
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', skewness=' num2str(skewness) ', runs=' num2str(runs)])
saveas(gcf,['plot/100_runs/binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% solid
load(['data/100_runs/binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '.mat'],'Linkmin','Cost', 'Load', 'K','K0'); 

Linkmin_Gbps = [7.5, 15, 22.5, 30];

figure;
face_alpha=1;
bar(Linkmin_Gbps,Cost(:,[3 1 2 4]),w1,'FaceAlpha',face_alpha);
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR')
xlabel('physical link capacity (Gbps)');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
saveas(gcf,['plot/100_runs/solid_binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_cost'],'epsc')

figure;
bar(Linkmin_Gbps,Load(:,[3 1 2 4 ]),w1,'FaceAlpha',face_alpha);
%hold on;
%plot([min(Linkmin) max(Linkmin)],[1 1],'r--','LineWidth',2);
%hold off;
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
xlabel('physical link capacity (Gbps)');
ylabel('congestion');
ylim([0.5 max(Load(:,2))+.1])
set(gca, 'FontSize', fontsize);
% title(['Topology ', topology, ', skewness=' num2str(skewness) ', runs=' num2str(runs)])
saveas(gcf,['plot/100_runs/solid_binary_cache_capacity_linkcapacity_C_' num2str(C) '_' topology '_skewness_' num2str(skewness) '_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)


