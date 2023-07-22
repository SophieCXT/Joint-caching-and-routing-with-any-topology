function [ cost, load, paths_output, links, cost_splittable, load_vector_splittable ] = bicriteria_MSUFP( G, c_link, s, d, lambda, K )
%Algorithm 2: Bicriteria Approximation for MSUFP
% Input:
% G: adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links).
% c_link: matrix containing link capacities; c_link(i,j) is the capacity of
% link (i,j). 
% s: source (single source)
% d: d(i) is destination for commodity i
% lambda: lambda(i) is demand for commodity i
% K: design parameter (assumed to be a positive integer)
% Output:
% cost: total routing cost
% load: E*1 array; load(e) load factor (flow/capacity) of link e
% paths_output: 1*n cell array; paths_output{i} is the path for commodity i as node sequence
% (from s to d(i))
%
% Note: need to run 'addpath ./matlab_bgl' before running this function (as it
% invokes dfs function from matlab_bgl). 

V = length(G); % #nodes
n = length(d); % #commodities
E = nnz(G); % #directed links
links = zeros(E,2); % links(j,:): the two endpoints of link j
link_index = zeros(size(G)); % link_index(u,v): index of link (u,v)
j = 0;
for u=1:V
    for v=1:V
        if G(u,v) ~= 0 
            j = j + 1;
            links(j,:) = [u v];
            link_index(u,v) = j;
        end
    end
end
        
%% line 1:
n_large = 10^6; % threshold for large matrix/vector
if max(E^2*n,V*E*n^2) <= n_large
    A = zeros(E,E*n);
    Aeq = zeros(V*n,E*n);    
else
    A = sparse(E,E*n);
    Aeq = sparse(V*n,E*n);
end
f = zeros(E*n,1);
for i=1:n
    for e=1:E
        f((i-1)*E+e) = lambda(i)*G(links(e,1),links(e,2));
    end
end
b = zeros(E,1);
for e=1:E
    A(e,([1:n]-1)*E+e) = reshape(lambda,1,n);
    b(e) = c_link(links(e,1),links(e,2));
end % A*X <= b (link capacity constraint)
for i=1:n
    for u=1:V
        for w=1:V
            if link_index(u,w) > 0
                Aeq((i-1)*V+u, (i-1)*E+link_index(u,w)) = 1;
            end
            if link_index(w,u) > 0
                Aeq((i-1)*V+u, (i-1)*E+link_index(w,u)) = -1;
            end
        end
    end
end
beq = zeros(V*n,1);
beq(([1:n]-1).*V+s) = 1;
for i=1:n
    beq((i-1)*V+d(i)) = -1;
end% Aeq*X = beq (flow conservation constraint)
lb = zeros(E*n,1);
ub = ones(E*n,1);
options = optimoptions('linprog','Display','none');
[X,fval] = linprog(f,A,b,Aeq,beq,lb,ub,options);
cost_splittable = fval; 
var_f = reshape(X,E,n); % var_f(e,i): fraction of commodity i routed on link e

f_splittable = var_f*lambda; % absolute traffic load on each link as E*1 vector
load_vector_splittable = f_splittable./b; 

%% line 2:
var_f_path = zeros(E,n); % var_f_path(p,i): fraction of commodity i routed on path p
paths = cell(E,n); % paths{p,i}: p-th path for commodity i as node sequence
costs = zeros(E,n); % costs(p,i): cost of paths{p,i}
for i=1:n
    A = zeros(V);
    for e=1:E
        A(links(e,1),links(e,2)) = var_f(e,i);
    end
    [dist, dt, ft, pre] = dfs(sparse(A),s); 
    p = 0;
    while dt(d(i)) ~= -1% find a new path
        p = p + 1;
        paths{p,i} = [pre(d(i)) d(i)];
        costs(p,i) = G(pre(d(i)), d(i)); 
        var_f_path(p,i) = A(pre(d(i)),d(i));
        v = d(i);
        while pre(v) ~= s
            v = pre(v);
            paths{p,i} = [pre(v) paths{p,i}];
            costs(p,i) = costs(p,i) + G(pre(v),v);
            var_f_path(p,i) = min(var_f_path(p,i),A(pre(v),v));
        end% paths{p,i} a new s->d(i) path, var_f_path(p,i) the "bottleneck capacity" on this path
        for j=1:length(paths{p,i})-1
            A(paths{p,i}(j),paths{p,i}(j+1)) = A(paths{p,i}(j),paths{p,i}(j+1)) - var_f_path(p,i); 
        end% update residual capacity
        [dist, dt, ft, pre] = dfs(sparse(A),s);
    end          
end
%% line 3:
lambda_bar = zeros(size(lambda));
lambda_max = max(lambda);
lambda_bar(lambda<lambda_max) = lambda_max.*2.^(floor(K.*log2(lambda(lambda<lambda_max)./lambda_max))./K);
lambda_bar(lambda==lambda_max) = lambda_max/2^(1/K); 
for i=1:n
    [~,I] = sort(costs(:,i),'descend');
    delta = 0; % amount of reduction in flow for commodity i
    j = 1;
    while delta < 1-lambda_bar(i)/lambda(i)
        delta1 = min(1-lambda_bar(i)/lambda(i)-delta, var_f_path(I(j),i));
        var_f_path(I(j),i) = var_f_path(I(j),i) - delta1; % reduce flow in descending order of path costs to satisfy the reduced demand lambda_bar
        delta = delta + delta1; 
        j = j + 1;
    end    
end
%% line 4:
S = false(K,n); % S(j,:): logicals of whether each commodity belongs to S_{j-1}
for j=1:K
    S(j,:) = (round(K.*(log2(lambda_max./lambda_bar)-floor(log2(lambda_max./lambda_bar)))) == K-j+1 | ...
        round(K.*(log2(lambda_max./lambda_bar)-floor(log2(lambda_max./lambda_bar)))) == -j+1);
end
% sanity check: every commodity belongs to one and only one S(j,:)
for i=1:n
    if sum(S(:,i)) ~= 1
        error(['commodity ' num2str(i) ' not included in any S_j']);
    end
end

%% line 5:
% var_f_path(:,S(j,:)) is the reduced flow for commodities in S_{j-1}, the
% corresponding paths are stored in paths{:,S(j,:)}. 
%% line 6-7:
paths_output = cell(1,n);
for j=1:K
    var_f_link = zeros(V); % link-level flow for commodities in S_{j-1}; var_f_link(u,v): (sum) flow (in absolute rate) on link (u,v) for var_f_path(:,S(j,:)) (in fraction of commodity)
    for i=1:n
        if S(j,i)
            for p=1:E
                for k=1:length(paths{p,i})-1
                    var_f_link(paths{p,i}(k),paths{p,i}(k+1)) = var_f_link(paths{p,i}(k),paths{p,i}(k+1)) + lambda(i)*var_f_path(p,i);
                end
            end
        end
    end
    % disp(['j = ' num2str(j) ', calling rounding_MSUFP...'])
    [ paths_output(S(j,:)) ] = rounding_MSUFP( G, s, d(S(j,:)), lambda_bar(S(j,:)), var_f_link, links, link_index  );
end

%% performance check:
cost = 0; % total routing cost under the computed solution
total_flow = zeros(E,1);
c_link_vector = zeros(E,1);
for i=1:n
    for k=1:length(paths_output{i})-1
        total_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))) = total_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))) + lambda(i);
        cost = cost + lambda(i)*G(paths_output{i}(k),paths_output{i}(k+1));    
    end
end
for e=1:E
    c_link_vector(e) = c_link(links(e,1),links(e,2));
    if total_flow(e) >= lambda_max*2^(1/K)/2/(2^(1/K)-1) + 2^(1/K)*c_link_vector(e)
        error(['link ' num2str(e) ' has load ' num2str(total_flow(e)) ', violating strict upper bound ' num2str(lambda_max*2^(1/K)/2/(2^(1/K)-1) + 2^(1/K)*c_link_vector(e))]);
    end
end
load = total_flow./c_link_vector; % load factor per link
disp(['MSUFP:: total cost: ' num2str(cost) ' (minimum cost of splittable flow: ' num2str(cost_splittable) '), maximum load factor: ' num2str(max(load)) ' (maximum load of splittable flow: ' num2str(max(load_vector_splittable)) ')' ])


end

