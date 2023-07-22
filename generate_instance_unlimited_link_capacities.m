function [  G,server,client,request ] = generate_instance_unlimited_link_capacities( A,cmin,cmax,cmin_s,cmax_s,deg_client,n_client,C,skewness )
%% generate a random instance of the problem (link costs, request rates):
V = length(A); 
% generate link costs:
G = zeros(size(A));
G(A>0) = cmin + (cmax-cmin).*rand(nnz(A),1);
% select client/server:
degree = sum(A,2); 
[~,I] = sort(degree); 
server = I(1); % server is the first node with the minimum degree %find(degree<=deg_server,1,'first'); 
% client = setdiff(find(degree<=deg_client),server);
client = I(2:n_client+1);
G(server,G(server,:)>0) = cmin_s + (cmax_s-cmin_s).*rand(1,sum(G(server,:)>0)); % server is remote, connected to test of network by links of delays in [100, 200]ms
% generate demands:
request_method = 3;
switch request_method
    case 1
        % method 1: Zipf over all C*#client request types
        request1 = [1:C*length(client)].^(-skewness);
        request1 = length(client).* request1 ./ sum(request1);
        request1 = request1(randperm(C*length(client))); % randomly permuted lambda_{i,s} with total rate = #clients
        request = zeros(C,V);
        request(:,client) = reshape(request1,C,length(client));
    case 2
        % method 2: Zipf over the C contents, equal for all clients
        request1 = [1:C].^(-skewness);
        request1 = length(client).* request1 ./ sum(request1); % request1(i): request rate for content i over all the clients
        request = zeros(C,V);
        for i=1:C
            request(i,client) = (request1(i)/length(client))*ones(1,length(client));
        end
    case 3
        % method 3: Zipf over the C contents at each client, with randomly
        % permuted ranking for each client
        request1 = [1:C].^(-skewness);
        request1 = request1' ./ sum(request1); % rates for requesting contents 1,...,C at each client (sum to one)
        request = zeros(C,V);
        for s=1:length(client)
            request(:,client(s)) = request1(randperm(C));
        end  
end

if false
%% plot the topology to find out where to put the server and the clients:
smin = 4; smax = 10;
g = graph(A);
figure;
h = plot(g, 'Layout', 'force');
highlight(h,server,'NodeColor','g','MarkerSize',smin+(smax-smin).*1); % green for server
highlight(h,client,'NodeColor','r','MarkerSize',smin+(smax-smin).*.5); % red for clients
switches = setdiff(1:V,[server; client]);
highlight(h,switches,'NodeColor','b', 'MarkerSize',smin); % blue for switches (node size indicates cache size)
end


end

