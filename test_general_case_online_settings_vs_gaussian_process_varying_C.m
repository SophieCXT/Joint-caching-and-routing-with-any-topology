% General case: finite cache capacities, finite link capacities
topology = 'Abovenet';
% topology = 'TiscaliEurope'; % not used
% topology = 'Level3';
% topology = 'VerioUS';

% addpath("matlab_bgl");
% addpath("C:\gurobi912\win64\matlab\");

switch topology
    case 'Abovenet'
        C_list = [10 20 30 40];
        C = 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]
%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]
        linkmin = 13715.796; % 5=n_client,minimum link capacity (before augmentation)
        linkmax = 13715.796; % maximum link capacity (before augmentation)
        Linkmin = [2 3 4 5].*13715.796./3; %5=n_client,
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
starting_point = 470; %start testing at 469+1, make sure we assumes the computation cost delay within 1 hours; New section: ‘Data driven simulation’,
weeks = 551;%hourly things ;weeks: 18;
predict_batch_size = 5; % one hour = one Monte Carlo run, predicting 5 hours at a time
predict_batch_num = 20; % num of iterations, GP_top30_videos_model_plus_5_pred_11_29_[1~20]

w1 = .85; 
w2 = .4;
total_num_request=34289494;

%% vary C :
iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(Linkmin),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
batch_num=0;
n_timestamp=0;
for i_p = 1:length(C_list)
    C = C_list(i_p);
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_online_settings( n_timestamp,batch_num, A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths);

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
Cost_online_settings = Cost./(predict_batch_num*predict_batch_size*total_num_request*1000*3600);
Load_online_settings = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_C_list_' topology '_linkmin_' num2str(linkmin) '_online_settings.mat'],'C_list','Cost_online_settings', 'Load_online_settings', 'iter_maximum', 'c_client'); 

% load(['data/general_case_linkcapacity_' topology '.mat']);
% figure;
% bar(Linkmin,Cost_online_settings);
% legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
% xlabel('catalog size');
% ylabel('routing cost');
% set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_online_settings_cost'],'epsc')
% figure;
% bar(Linkmin,Load_online_settings);
% legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
% xlabel('catalog size');
% ylabel('congestion');
% % title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
% set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_online_settings_congestion'],'epsc')
% 
% % add voice when stop
% load handel;
% player = audioplayer(y, Fs);
% play(player)
%% gaussian_process: vary C:

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(Linkmin),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
batch_num=0;
n_timestamp=0;
for i_p = 1:length(C_list)
    C = C_list(i_p);
for batch_num=1:predict_batch_num
for n_timestamp=1:predict_batch_size
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general_gaussian_process( n_timestamp,batch_num, A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths);
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
Cost_gaussian_process = Cost./(predict_batch_num*predict_batch_size*total_num_request*1000*3600);%runs;
Load_gaussian_process = Load./(predict_batch_num*predict_batch_size);
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs_trace_driven/general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_gaussian_process.mat'],'C_list','Cost_gaussian_process', 'Load_gaussian_process', 'iter_maximum', 'c_client'); 

% load(['data/general_case_linkcapacity_' topology '.mat']);
% figure;
% bar(Linkmin,Cost_gaussian_process);
% legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
% xlabel('catalog size');
% ylabel('routing cost');
% set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_cost'],'epsc')
% figure;
% bar(Linkmin,Load_gaussian_process);
% legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
% xlabel('catalog size');
% ylabel('congestion');
% % title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
% set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_gaussian_process_congestion'],'epsc')
% 
% % add voice when stop
% load handel;
% player = audioplayer(y, Fs);
% play(player)


%% plot online_vs_gaussian_process_
% vary C_list 
% running sequence: first part, then gaussian_process vary (minimum) link capacity,
% then first part, then 'vary (minimum) link capacity'
figure;
b=bar(C_list,Cost_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(C_list,Cost_gaussian_process,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('catalog size');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_cost'],'epsc')
% saveas(gcf,['plot/predict_window_5/online_vs_gaussian_process_general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_cost'],'fig')


figure;
b=bar(C_list,Load_online_settings,w1,'FaceAlpha',0.5);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
hold on
%gaussian_process
b=bar(C_list,Load_gaussian_process,w2,'FaceAlpha',0.9);
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4940, 0.1840, 0.5560];%purple
b(5).FaceColor = [0.4660, 0.6740, 0.1880]; %green, 
grid on
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('catalog size');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
% saveas(gcf,['plot/100_runs_trace_driven/online_vs_gaussian_process_general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_congestion'],'epsc')
% saveas(gcf,['plot/predict_window_5/online_vs_gaussian_process_general_case_C_list_linkmin_' num2str(linkmin) '_' topology '_cost'],'fig')

