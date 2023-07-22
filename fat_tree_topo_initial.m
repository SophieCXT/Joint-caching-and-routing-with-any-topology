function [G,c_link] = fat_tree_topo_initial(N,Pe,Pc,Bl)
% fat_tree_topo_initial is a function initializing the topology 2 - level fat-tree, calculating
% the number of nodes, core switch number, edge switch number, link, etc.
% And the input contains:
%  eg, N=128 nodes, with Pe=Pc=36 ports and a blocking factor Bl=1   
%  Eptn=Trunc (Pe * Bl / (1+Bl))
%  Eptc=Pe – Eptn
%  using Bl=1 for non-blocking networks
Eptn=floor (Pe * Bl / (1+Bl)); % Eptn — “edge ports to nodes”
Eptc=Pe-Eptn; % Eptc — “edge ports to core”:
E=ceil(N/Eptn); % E: number of edge switch
B=fixDiv(Pc,E);
C=ceil(Eptc/B);

% total number of nodes: N + E + C
G = zeros(N + E + C); % create a graph that contains N + E + C nodes: host/server nodes, edge switch nodes, core switch nodes
% for load balance, each edge switch will connect ceil(N/E) hosts
num_host = ceil(N/E);
% for the first level
for i=1:E
    if i == E
        G(num_host*(i-1)+1:N,N+i)=1; % set virtual link costs = 1 (as G(.,.)=0 is interpreted as no link)
        % the final cost needs to be reduced by sum(sum(request))*1 to ignore this
        % cost of virtual links. 
    else
        G(num_host*(i-1)+1:num_host*i,N+i)=1;
    end  
end
%for the second level: 
num_edgeswitch = ceil(E/C);
for i=1:C
    if i == C
        G(N+num_edgeswitch*(i-1)+1:E,N++Ei)=1; % set virtual link costs = 1 (as G(.,.)=0 is interpreted as no link)
        % the final cost needs to be reduced by sum(sum(request))*1 to ignore this
        % cost of virtual links. 
    else
        G(N+num_edgeswitch*(i-1)+1:num_edgeswitch*i,N+E+i)=1;
    end  
end

G = max(G,G'); %symmetric attributes for a graph,bidirection link 
c_link = zeros(size(G));
%link cost initialization
end

