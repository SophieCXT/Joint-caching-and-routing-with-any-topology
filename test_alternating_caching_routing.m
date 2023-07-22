% toy example for debugging "source_selection_routing_given_placement" and "integral_caching_given_paths":
c_infty = 1000;
c_epsilon = .1;
c_large = 10;
b_epsilon = 1;
b_large = 100;
C = 2; 
G = zeros(6);
c_link = zeros(6);
c_link(1,4) = C*b_epsilon; G(1,4) = c_epsilon; 
c_link(2,4) = b_epsilon; G(2,4) = c_epsilon;
c_link(3,4) = b_epsilon; G(3,4) = c_large;
c_link(4,5) = b_epsilon; G(4,5) = c_epsilon;
c_link(3,5) = b_epsilon; G(3,5) = c_epsilon;
c_link(4,6) = b_large; G(4,6) = c_infty; 
G = max(G,G');
c_link = max(c_link,c_link');
c_cache=[0 2 1 0 0 2]';
request=zeros(C,6);
request(1,1)=.9;
request(2,1)=.8;

placement = zeros(6,C); placement(6,:)=1; % initially only placed at the remote server
integral = 1; % integral routing (i.e., single-path routing)?
heuristic = 1; % compute routing by heuristic (only if integral routing)?
cost_final = inf; 
tic
for i=1:5
    disp(['iteration ' num2str(i) ':'])
    [ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement, request, integral, heuristic );
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing routing'])
    [ cost, placement ] = integral_caching_given_paths( G, c_cache, request, flows, paths, request_type ); 
    disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after optimizing caching'])
    disp(' ')
    if cost >= cost_final
        disp(['no improvement after ' num2str(i) ' iterations; stop'])
        break;
    else
        cost_final = cost;
    end
end
disp(['alternating optimization takes ' num2str(toc) ' sec'])
%%
placement = zeros(6,C); placement(6,:)=1; % initially only placed at the remote server
tic
[ cost, load_max, placement, flows, paths, request_type  ] = greedy_caching( G, c_cache, c_link, request, integral, heuristic, placement );
disp(['cost = ' num2str(cost) ', congestion = ' num2str(load_max) ' after greedy caching'])
disp(['greedy optimization takes ' num2str(toc) ' sec'])
% Observation: alternating optimization is faster than greedy content
% placmement (based on optimal/heuristic routing)

c_opt = request(1,1)*(G(2,4)+G(4,1))+request(2,1)*(G(3,5)+G(5,4)+G(4,1)); 
var_x_opt = [0 0; 1 0; 0 1; 0 0; 0 0; 1 1];
