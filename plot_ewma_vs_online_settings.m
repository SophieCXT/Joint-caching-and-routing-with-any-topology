% Overlay Bar Graphs
%% case one: unlimited link capacities
% Special case 1 of online settings: unlimited link capacities (routing = RNR, integral content placement)
% Data is from youtube viewing history
topology = 'Abovenet';
% topology = 'Level3';
% topology = 'VerioUS';

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
runs = 1; % tian, at least ten monte carlo runs.
window_size = 2;
weeks = 3;
fontsize = 12;
w1 = .85; 
w2 = .4;

% load data from online settings and ewma predictions

%% vary edge cache capacity:
% load(['data/online_settings/no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology '.mat']); 
% load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\no_link_capacity_candidatepaths_C_10_edgecachesize_2_Abovenet.mat')
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\no_link_capacity_cachecapacity_C_10_edgecachesize_4_Abovenet.mat')
figure;
% bar(C_client, C_RNR,w1,'facecolor',color_matrix(1:2:5,:));
b=bar(C_client, C_RNR,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
hold on
% saveas(gcf,['data/online_settings/no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')
%ewma
% load(['data/online_settings/ewma_no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology '.mat']); 
% load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\ewma_no_link_capacity_candidatepaths_C_10_edgecachesize_2_Abovenet.mat')
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\ewma_no_link_capacity_cachecapacity_C_10_edgecachesize_4_Abovenet.mat')
% bar(C_client, C_RNR,w1,'facecolor',color_matrix(2:2:6,:));
b=bar(C_client, C_RNR,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
legend('Alg. 1','k shortest paths','shortest path','EWMA-Alg. 1','EWMA-k shortest paths','EWMA-shortest path');
grid on
ylabel('routing cost');
xlabel('cache capacity');
set(gca, 'FontSize', fontsize);
saveas(gcf,['data/online_settings/ewma_vs_online_no_link_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% vary k_paths:
% load(['data/online_settings/no_link_capacity_candidatepaths_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology '.mat'], 'K_paths','C_RNR','c_client'); 
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\no_link_capacity_candidatepaths_C_10_edgecachesize_2_Abovenet.mat')
figure;
b=bar(K_paths, C_RNR,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
% b(3).FaceColor = [0.9290 0.6940 0.1250];
hold on
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\ewma_no_link_capacity_candidatepaths_C_10_edgecachesize_2_Abovenet.mat')
b=bar(K_paths, C_RNR,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
% b(3).FaceColor = [0.9290 0.6940 0.1250];
legend('Alg. 1','k shortest paths','EWMA-Alg. 1','EWMA-k shortest paths')
xlabel('#candidate paths');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
saveas(gcf,['data/online_settings/ewma_vs_online_no_link_capacity_candidatepaths_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')

% % add voice when stop
% load handel;
% player = audioplayer(y, Fs);
% play(player)

%% plot ewma vs online settings
Linkmin = [.2 .3 .4 .5]; %[.5 .75 1 1.25 ];
K_total = [2:2:50]; % optimal K = 12
K = 32; %12; % optimized for ours 

fontsize = 12;
w1 = .85; 
w2 = .4;
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\binary_cache_capacity_K_C_10_Abovenet.mat') 
figure;
% bar(C_client, C_RNR,w1,'facecolor',color_matrix(1:2:5,:));
b=bar(Linkmin,Load(:,[3 1 2 4 ]),w1);%,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
hold on
%ewma
load('C:\Users\tbx5027\The Pennsylvania State University\He, Ting - SDN security\Ting team\joint placement, selection, routing\matlab - Tian\data\online_settings\binary_cache_capacity_linkcapacity_C_10_Abovenetewma.mat') %binary_cache_capacity_K_C_10_Abovenet
bar(Linkmin,Load(:,[3 1 2 4 ]),w2);
% b=bar(Linkmin,Load(:,[3 1 2 4 ]))%;,w2,'FaceAlpha',0.9);%bar(C_client, C_RNR,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660, 0.6740, 0.1880]; %green
legend('splittable flow', ['Alg. 2: K = ' num2str(K)],['Alg. 2: K = ' num2str(K0)], 'RNR'); %, ['bound: K = ' num2str(K)],['bound: K = ' num2str(K0)])
grid on
xlabel('link capacity');
ylabel('congestion');
% ylim([0.5 max(Load(:,2))+.1])
set(gca, 'FontSize', fontsize);
saveas(gcf,['data/online_settings/ewma_vs_online_binary_capacity_cachecapacity_C_' num2str(C) '_edgecachesize_' num2str(c_client) '_' topology],'epsc')
