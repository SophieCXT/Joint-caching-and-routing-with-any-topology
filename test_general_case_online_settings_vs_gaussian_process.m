% General case: finite cache capacities, finite link capacities
topology = 'Abovenet';
% topology = 'TiscaliEurope'; % not used
% topology = 'Level3';
% topology = 'VerioUS';

% addpath("matlab_bgl");
% addpath("C:\gurobi912\win64\matlab\");
% >> cd C:\gurobi912\win64\matlab
% >> gurobi_setup

switch topology
    case 'Abovenet'
        C = 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; 
        linkmin = 13715.796; % .5; % minimum link capacity (before augmentation)
        linkmax = 13715.796; %equals 15 Gbps % maximum link capacity (before augmentation)
        Linkmin = [2 3 4 5].*13715.796./3;%5=n_client, [.1 .2 .3 .4] [.2 .3 .4 .5]; %[.5 .75 1 1.25 ];
        Linkmin_Gbps = [10 15 20 25]; %in Gbps
       
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;
        deg_client = 3; % maximum degree of clients
        n_client = 5;
        k_paths = 10; % [Ioannidis18JSAC]
    case 'TiscaliEurope'
        C = 10;
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]
        linkmin = .2; % minimum link capacity (before augmentation)
        linkmax = .2; % maximum link capacity (before augmentation)
        Linkmin = [.2 .3 .4 .5];
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;
        deg_client = 1; % maximum degree of clients
        n_client = 12;
        k_paths = 30; % [Ioannidis18JSAC]
    case 'Level3'
        C = 15; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        n_client = 8;
        linkmin = .2./n_client; %.4; % minimum link capacity (before augmentation)
        linkmax = .2./n_client; %.4; % maximum link capacity (before augmentation)
        Linkmin = [.2 .4 .6 .8]./n_client; %[.5 .75 1 1.25 ];
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;        
        deg_client = 5; % maximum degree of clients
        k_paths = 10; % [Ioannidis18JSAC]
    case 'VerioUS'
        C = 30; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)        
        skewness = 1.2; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]        
        linkmin = .2; % minimum link capacity (before augmentation)
        linkmax = .2; % maximum link capacity (before augmentation)
        Linkmin = [.2 .4 .6 .8]; %[.5 .75 1 1.25 ];
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;                 
        deg_client = 1; %17; % maximum degree of clients
        n_client = 12; % #clients (with degree up to deg_client); note: this includes caches (but excludes the server)
        k_paths = 10; % [Ioannidis18JSAC]
end

max_iter = 10; % max #iterations (to bound running time) for alternating optimization
load(['data/' topology '.mat']); % A: sparse adjacency matrix with binary entries
A = full(A);
V = length(A); 
runs = 10;
fontsize = 16; %used to be 12
window_size = 2; %for gaussian_process
weeks = 551;%hourly things ;weeks: 18;
predict_batch_size = 5; % one hour = one Monte Carlo run, predicting 5 hours at a time
predict_batch_num = 20; % num of iterations, GP_top30_videos_model_plus_5_pred_11_29_[1~20]

w1 = .85; 
w2 = .4;

total_num_request=34289494;


%% vary link capacity:
% Linkmin = [.2 .4 .6 .8];
alpha = 1; % linkmax/linkmin = alpha --> all links have identical capacities before augmentation

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(Linkmin),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_online_settings( n_timestamp,batch_num, A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths);
for i_p = 1:length(Linkmin)
    linkmin = Linkmin(i_p);
    linkmax = alpha*linkmin;
    disp(' ')
    disp(['link capacity (before augmentation) = ' num2str(linkmin) ':'])
    disp(' ')
    % overwrite link capacities, but retain other parameters:
    c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link;
    c_link = max(c_link,demand_SP);
[ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths );
Cost(i_p,1) = Cost(i_p,1)+cost_alternating_FR;
Load(i_p,1) = Load(i_p,1)+load_alternating_FR;
Cost(i_p,2) = Cost(i_p,2)+cost_alternating; 
Load(i_p,2) = Load(i_p,2)+load_alternating; 
Cost(i_p,3) = Cost(i_p,3)+cost_SP;
Load(i_p,3) = Load(i_p,3)+load_SP;
Cost(i_p,4) = Cost(i_p,4)+cost_SP2;
Load(i_p,4) = Load(i_p,4)+load_SP2; 
Cost(i_p,5) = Cost(i_p,5)+cost_kSP_RNR;
Load(i_p,5) = Load(i_p,5)+load_kSP_RNR;
iter_maximum = max([iter_maximum, iter_alternating, iter_alternating_FR]);
end
end
end
Cost_online_settings = Cost./total_num_request;
Load_online_settings = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C) '_' topology '_linkmin_' num2str(linkmin) '_online_settings.mat'],'Linkmin','Cost', 'Load', 'iter_maximum', 'c_client'); 

% load(['data/general_case_linkcapacity_' topology '.mat']);
figure;
bar(Linkmin_Gbps,Cost_online_settings);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('link capacity');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_online_settings_cost'],'epsc')
figure;
bar(Linkmin_Gbps,Load_online_settings);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('link capacity');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C)  '_linkmin_' num2str(linkmin) '_' topology '_online_settings_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)
%% gaussian_process: vary link capacity:
% Linkmin = [.2 .4 .6 .8];
alpha = 1; % linkmax/linkmin = alpha --> all links have identical capacities before augmentation

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(Linkmin),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
% 
% n_timestamp=0;
% batch_num=0;

for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_gaussian_process( n_timestamp,batch_num, A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths);
for i_p = 1:length(Linkmin)
    linkmin = Linkmin(i_p);
    linkmax = alpha*linkmin;
    disp(' ')
    disp(['link capacity (before augmentation) = ' num2str(linkmin) ':'])
    disp(' ')
    % overwrite link capacities, but retain other parameters:
    c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link;
    c_link = max(c_link,demand_SP);
[ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths );
Cost(i_p,1) = Cost(i_p,1)+cost_alternating_FR;
Load(i_p,1) = Load(i_p,1)+load_alternating_FR;
Cost(i_p,2) = Cost(i_p,2)+cost_alternating; 
Load(i_p,2) = Load(i_p,2)+load_alternating; 
Cost(i_p,3) = Cost(i_p,3)+cost_SP;
Load(i_p,3) = Load(i_p,3)+load_SP;
Cost(i_p,4) = Cost(i_p,4)+cost_SP2;
Load(i_p,4) = Load(i_p,4)+load_SP2; 
Cost(i_p,5) = Cost(i_p,5)+cost_kSP_RNR;
Load(i_p,5) = Load(i_p,5)+load_kSP_RNR;
iter_maximum = max([iter_maximum, iter_alternating, iter_alternating_FR]);
end
end
end
Cost_gaussian_process = Cost./total_num_request;%runs;
Load_gaussian_process = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process.mat'],'Linkmin','Cost_gaussian_process', 'Load_gaussian_process', 'iter_maximum', 'c_client'); 

% load(['data/general_case_linkcapacity_' topology '.mat']);
figure;
bar(Linkmin_Gbps,Cost_gaussian_process);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('link capacity');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_cost'],'epsc')
figure;
bar(Linkmin_Gbps,Load_gaussian_process);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('link capacity');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs_trace_driven/general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% vary cache capacity:
% C_client = 1:4; 

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(C_client),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_online_settings( n_timestamp,batch_num, A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, min(C_client), k_paths);
for i_p = 1:length(C_client)
    c_client = C_client(i_p);
    disp(['c_client = ' num2str(c_client) ':'])
    disp(' ')
    % overwrite cache capacities, but retain other parameters:
    c_cache(client) = c_client;
[ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths );
Cost(i_p,1) = Cost(i_p,1)+cost_alternating_FR;
Load(i_p,1) = Load(i_p,1)+load_alternating_FR;
Cost(i_p,2) = Cost(i_p,2)+cost_alternating; 
Load(i_p,2) = Load(i_p,2)+load_alternating; 
Cost(i_p,3) = Cost(i_p,3)+cost_SP;
Load(i_p,3) = Load(i_p,3)+load_SP;
Cost(i_p,4) = Cost(i_p,4)+cost_SP2;
Load(i_p,4) = Load(i_p,4)+load_SP2; 
Cost(i_p,5) = Cost(i_p,5)+cost_kSP_RNR;
Load(i_p,5) = Load(i_p,5)+load_kSP_RNR;
iter_maximum = max([iter_maximum, iter_alternating, iter_alternating_FR]);
end
end
end
Cost_online_settings = Cost./total_num_request;
Load_online_settings = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '.mat'],'C_client','Cost_online_settings', 'Load_online_settings', 'iter_maximum', 'linkmin','linkmax'); 

figure;
bar(C_client,Cost_online_settings);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
ylabel('routing cost/request (ms)');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_cost'],'epsc')
figure;
bar(C_client,Load_online_settings);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
ylabel('routing cost/request (ms)');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% gaussian_process: vary cache capacity:
% C_client = 1:4; 

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(C_client),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 

n_timestamp=0;
batch_num=0;

for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_gaussian_process( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, min(C_client), k_paths);
for i_p = 1:length(C_client)
    c_client = C_client(i_p);
    disp(['c_client = ' num2str(c_client) ':'])
    disp(' ')
    % overwrite cache capacities, but retain other parameters:
    c_cache(client) = c_client;
[ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths );
Cost(i_p,1) = Cost(i_p,1)+cost_alternating_FR;
Load(i_p,1) = Load(i_p,1)+load_alternating_FR;
Cost(i_p,2) = Cost(i_p,2)+cost_alternating; 
Load(i_p,2) = Load(i_p,2)+load_alternating; 
Cost(i_p,3) = Cost(i_p,3)+cost_SP;
Load(i_p,3) = Load(i_p,3)+load_SP;
Cost(i_p,4) = Cost(i_p,4)+cost_SP2;
Load(i_p,4) = Load(i_p,4)+load_SP2; 
Cost(i_p,5) = Cost(i_p,5)+cost_kSP_RNR;
Load(i_p,5) = Load(i_p,5)+load_kSP_RNR;
iter_maximum = max([iter_maximum, iter_alternating, iter_alternating_FR]);
end
end
end
Cost_gaussian_process = Cost./total_num_request;
Load_gaussian_process = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process.mat'],'C_client','Cost_gaussian_process', 'Load_gaussian_process', 'iter_maximum', 'linkmin','linkmax'); 

figure;
bar(C_client,Cost_gaussian_process);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
ylabel('routing cost/request (ms)');
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_cost'],'epsc')
figure;
bar(C_client,Load_gaussian_process);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
ylabel('routing cost/request (ms)');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_congestion'],'epsc')

% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)

%% plot online_vs_gaussian_process_
% vary link capacity 
% running sequence: first part, then gaussian_process vary (minimum) link capacity,
% then first part, then 'vary (minimum) link capacity'
figure;
b=bar(Linkmin_Gbps,Cost_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(Linkmin_Gbps,Cost_gaussian_process,w2,'FaceAlpha',1);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('physical link capacity (Gbps)');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_cost'],'epsc')


figure;
b=bar(Linkmin_Gbps,Load_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(Linkmin_Gbps,Load_gaussian_process,w2,'FaceAlpha',1);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('physical link capacity (Gbps)');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_linkcapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_congestion'],'epsc')

%% Plot: vary cache capacity 
% running sequence: first part, then gaussian_process vary (minimum) link capacity,
% then first part, then 'vary (minimum) link capacity'
% Cost_gaussian_process = Cost_gaussian_process./17;
% Load_gaussian_process = Load_gaussian_process./17;
figure;
b=bar(C_client,Cost_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(C_client,Cost_gaussian_process,w2,'FaceAlpha',1);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
ylabel('routing cost/request (ms)');
xlabel('cache capacity');
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_cost'],'epsc')


figure;
b=bar(C_client,Load_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(C_client,Load_gaussian_process,w2,'FaceAlpha',1);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('cache capacity');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_cachecapacity_C_' num2str(C) '_linkmin_' num2str(linkmin) '_' topology '_congestion'],'epsc')

