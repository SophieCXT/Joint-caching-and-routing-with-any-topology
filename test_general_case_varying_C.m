% General case: finite cache capacities, finite link capacities
topology = 'Abovenet';
% topology = 'TiscaliEurope'; % not used
% topology = 'Level3';
% topology = 'VerioUS';

% addpath('matlab_bgl');

switch topology
    case 'Abovenet'
        C_list = [10 20 30 40];
        C = 10; % total #items, following [Ioannidis18JSAC]5
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)
        skewness = 1.2; % skewness in demands, with average rate of 1 per client [Ioannidis18JSAC]
%         skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]
        linkmin = .2; % minimum link capacity (before augmentation)
        linkmax = .2; % maximum link capacity (before augmentation)
        Linkmin = [.1 .2 .3 .4];%[.2 .3 .4 .5];
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
        linkmin = .2; %.4; % minimum link capacity (before augmentation)
        linkmax = .2; %.4; % maximum link capacity (before augmentation)
        Linkmin = [.2 .4 .6 .8]; %[.5 .75 1 1.25 ];
        c_v = 0; % cache capacity at switches/backbone nodes
        c_client = 2; % cache capacity at clients/edge nodes, following [Ioannidis18JSAC]
        C_client = 1:4;        
        deg_client = 5; % maximum degree of clients
        n_client = 8;
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
fontsize = 16;
face_alpha=1;%0.8;
w1 = .85; 
w2 = .4;

%% vary vary C :
% Linkmin = [.2 .4 .6 .8];
alpha = 1; % linkmax/linkmin = alpha --> all links have identical capacities before augmentation

iter_maximum = 0; % record the max #iterations used by alternating optimization in all settings and all Monte Carlo runs
Cost = zeros(length(Linkmin),6); % 1: FR alternating, 2: IR alternating, 3: IR shortest path [Ioannidis'18TON], 4: IR shortest path+RNR (special case of [Ioannidis'18JSAC] with one candidate path per client), 5: IR k-shortest paths+RNR ([Ioannidis'18JSAC] with k_paths candidate paths per client)
Load = zeros(size(Cost)); 
for i_p = 1:length(C_list)
    C = C_list(i_p);
for rindx = 1:runs
    [ G,server,client,request, c_cache, c_link, rand_c_link,demand_SP ] = generate_instance_general( A,cmin,cmax,cmin_s,cmax_s,C,skewness,deg_client,n_client, linkmin, linkmax, c_v, c_client, k_paths);
[ cost_alternating_FR, load_alternating_FR, cost_alternating, load_alternating, cost_SP, load_SP, cost_SP2, load_SP2, cost_kSP_RNR, load_kSP_RNR, iter_alternating, iter_alternating_FR ] = simulate_general_case( G,server,client,request, c_cache, c_link, max_iter, k_paths );
Cost(i_p,1) = Cost(i_p,1)+cost_alternating_FR./sum(sum(request));
Load(i_p,1) = Load(i_p,1)+load_alternating_FR;
Cost(i_p,2) = Cost(i_p,2)+cost_alternating./sum(sum(request)); 
Load(i_p,2) = Load(i_p,2)+load_alternating; 
Cost(i_p,3) = Cost(i_p,2)+cost_SP./sum(sum(request));
Load(i_p,3) = Load(i_p,3)+load_SP;
Cost(i_p,4) = Cost(i_p,2)+cost_SP2./sum(sum(request));
Load(i_p,4) = Load(i_p,4)+load_SP2; 
Cost(i_p,5) = Cost(i_p,2)+cost_kSP_RNR./sum(sum(request));
Load(i_p,5) = Load(i_p,5)+load_kSP_RNR;
iter_maximum = max([iter_maximum, iter_alternating, iter_alternating_FR]);
end
end
Cost = Cost./runs;
Load = Load./runs;
disp(['maximum #iterations for alternating optimization (IR or FR): ' num2str(iter_maximum)])
% Observation: alternating optimization converges within 5 iterations in
% all the cases

save(['data/100_runs/general_case_C_list_' topology '.mat'],'Linkmin','Cost', 'Load', 'iter_maximum', 'c_client'); 

% load(['data/general_case_linkcapacity_' topology '.mat']);
figure;
bar(C_list,Cost,w1,'FaceAlpha',face_alpha);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('catalog size');
ylabel('routing cost/request (ms)');
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs/solid_general_case_C_list_' topology '_cost'],'epsc')

figure;
bar(C_list,Load,w1,'FaceAlpha',face_alpha);
legend('IC-FR (alternating)','IC-IR: alternating','        SP','        SP + RNR','        k-SP + RNR')
xlabel('catalog size');
ylabel('congestion');
% title(['Topology ', topology, ', C=' num2str(C), ', skewness=' num2str(skewness)])
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs/solid_general_case_C_list_' topology '_congestion'],'epsc')


% add voice when stop
load handel;
player = audioplayer(y, Fs);
play(player)