function [ G,server,cache,client,request, c_link, rand_c_link,demand_SP ] = generate_instance_binary_cache_capacities_online_settings( n_timestamp,batch_num,A,cmin,cmax,cmin_s,cmax_s,C,skewness,n_cache,deg_client,n_client, linkmin, linkmax, cache_pre)
% 
% %Import data
% t = csvread('data/weeks_views.csv',1,0);% = readtable('data/weeks_views.csv');
% 
% for video = 1:C %50,length(weeks_views), top 50 (C)
%     views_series = t(video,:).';
%     exp_ts(video,:) = tsmovavg(views_series,'e',window_size,1); %ts as timestamp
% %      triangular moving average moving average.
% end
% disp(['EWMA Predicted request has been generated, where window_size n=' num2str(window_size)])

m = readmatrix(['data/top_300_raw_diff_views.csv']);%batch_num, GP_top30_videos_model_plus_5_pred_11_29_[1~20]; ref: https://ww2.mathworks.cn/help/matlab/ref/readmatrix.html

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
if cache_pre == 0 % if not predetermine the caches:
%     cache = setdiff(find(degree<=degree(I(n_cache+1))), server); cache = cache(randperm(length(cache))); cache = cache(1:n_cache); % caches are randomly chosen from the next min-degree nodes
    cache = I(2:n_cache+1); % caches are the next n_cache min-degree node %setdiff(find(degree<=deg_cache), server);
else
    cache = cache_pre;
end
% client = setdiff(find(degree<=deg_client),[server; cache]);
client = I(n_cache+2:n_client+1);
G(server,G(server,:)>0) = cmin_s + (cmax_s-cmin_s).*rand(1,sum(G(server,:)>0)); % server is remote, connected to test of network by links of delays in [100, 200]ms
% generate demands over the top C videos at each client
% request1 = t(1:C,n_week); % top C videos, #views in n_week
request1 = m(550+5*(batch_num-1)+n_timestamp,2:C+1).'; % top C videos, #views in n_timestamp in batch_num th prediction file, each time predict 5 hours
request_client = rand(C,length(client));% randomly generates request distribution. The rows denote the requested contents, the columns denote clients.
for i=1:C
    % The logic is that each sum(row m )= request1(row m), which
    % satisfied randomness of: 1.clients request rate 2.clients behavior
    request_client(i,:) = (request_client(i,:) ./ sum(request_client(i,:))) .* request1(i);
end
request = zeros(C,V);
for s=1:length(client)
    request(:,client(s)) = request_client(:,s);
end

% generate link capacities:
c_link = zeros(size(A));
linkmin = max(linkmin, max(request1)); % make sure that minimum link capacity >= maximum demand as typically assumed in the literature of MSUFP
rand_c_link = rand(nnz(A),1); % store the random numbers drawn for computing link capacities (before augmentation) to increase linkmax proportionally
c_link(A>0) = linkmin + (linkmax-linkmin).*rand_c_link; %rand(nnz(A),1);
% augment link capacities to ensure it is feasible to serve each request
% from the server along the shortest path:
[~,SP] = Dijkstra_source(G, server);
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
smin = 4; 
smax = 10;
g = graph(A);
figure;
h = plot(g, 'Layout', 'force');
highlight(h,server,'NodeColor','g','MarkerSize',smin+(smax-smin)*1); % green for server
highlight(h,cache,'NodeColor','c','MarkerSize',smin+(smax-smin)*.5); % cyan for caches
highlight(h,client,'NodeColor','r','MarkerSize',smin+(smax-smin)*0); % red for clients
switches = setdiff(1:V,[server; cache; client]); 
highlight(h,switches,'NodeColor','b', 'MarkerSize',smin); % blue for switches      
for s=1:V
    for t=1:V
        if A(s,t)>0
            highlight(h,s,t,'LineWidth', .5+ 2*c_link(s,t)/length(client));
        end
        if demand_SP(s,t)>0
            highlight(h,s,t,'EdgeColor','r');
        end
    end
end    

end


end

