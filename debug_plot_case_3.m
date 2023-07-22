w1 = .85; 
w2 = .4;
fontsize = 16;
total_num_request=34289494;
load('data/100_runs_trace_driven/general_case_linkcapacity_C_10_Abovenet_linkmin_22859.66_online_settings.mat')
Cost_online_settings = Cost./total_num_request;
Load_online_settings = Load./100;
load('data/100_runs_trace_driven/general_case_C_list_linkmin_13715.796_Abovenet_gaussian_process.mat')
Cost_gaussian_process = Cost_gaussian_process.*360000000;%runs;
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
ylabel('routing cost');
set(gca, 'FontSize', fontsize);
saveas(gcf,['plot/100_runs_trace_driven/solid_case3_vary_link_capacity_cost'],'epsc')


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
saveas(gcf,['plot/100_runs_trace_driven/solid_case3_vary_link_capacity_congestion'],'epsc')


load('data/100_runs_trace_driven/general_case_cachecapacity_C_10_linkmin_13715.796_Abovenet.mat')
Cost_online_settings=Cost_gaussian_process.*360000000;
load('data/100_runs_trace_driven/general_case_cachecapacity_C_10_linkmin_13715.796_Abovenet_gaussian_process.mat')
Cost_gaussian_process = Cost_gaussian_process.*360000000;%runs;

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
saveas(gcf,['plot/100_runs_trace_driven/solid_case3_vary_cache_capacity_cost'],'epsc')


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
saveas(gcf,['plot/100_runs_trace_driven/solid_case3_vary_cache_capacity_congestion'],'epsc')

