function [  G,server,client,request ] = online_settings_generate_instance_unlimited_link_capacities( n_week,window_size, A,cmin,cmax,cmin_s,cmax_s,n_client,C )

%Import data
t = csvread('data/weeks_views.csv',1,0);% = readtable('data/weeks_views.csv');

for video = 1:50%length(weeks_views), top 50 (C)
    views_series = t(video,:).';
    exp_ts(video,:) = tsmovavg(views_series,'e',window_size,1); %ts as timestamp
%      triangular moving average moving average.
end

disp(['EWMA Predicted request has been generated, where window_size n=' num2str(window_size)])

%% generate a random instance of the online settings problem according to youtube data (link costs, request rates):
V = length(A); 
% generate link costs: Fixed as constant 
G = zeros(size(A));
G(A>0) = cmin + (cmax-cmin).*rand(nnz(A),1);
% select client/server:
degree = sum(A,2); 
[~,I] = sort(degree); 
server = I(1); % server is the first node with the minimum degree %find(degree<=deg_server,1,'first'); 
% client = setdiff(find(degree<=deg_client),server);
client = I(2:n_client+1);
G(server,G(server,:)>0) = cmin_s + (cmax_s-cmin_s).*rand(1,sum(G(server,:)>0)); % server is remote, connected to test of network by links of delays in [100, 200]ms
% generate demands over the top C videos at each client
request1 = t(1:C,n_week); % top C videos, #views in n th week
request1 = request1 ./ sum(request1); % rates for requesting contents 1,...,C at each client (sum to one)
request_client = rand(C,length(client));% randomly generates request distribution. The rows denote the requested contents, the columns denote clients.
for i=1:C
    % The logic is that each sum(row m )= request1(row m), which
    % satisfied randomness of: 1.clients request rate 2.clients behavior
    request_client(i,:) = (request_client(i,:) ./ sum(request_client(i,:))) .* request1(i);
end
% If use uniform distribution
%  generate N random numbers in the interval (a,b) with the formula,
% r = a + (b-a).*rand(N,1), with mean (b+a)/2 and variance ((b-a)^2)/12. 

request = zeros(C,V);
for s=1:length(client)
    request(:,client(s)) = request_client(:,s);
end



