function [ cost, load_max, flows, paths, request_type ] = source_selection_routing_given_placement( G, c_link, placement, request, integral, heuristic )
%Source selection and routing under a given integral content placement:
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_link: V*V matrix containing link capacities; c_link(i,j) is the capacity of
% link (i,j). 
% placement: V*C matrix; placement(v,i) = x_{vi} (whether node v stores
% item i). 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% integral: 1 if integral source selection & routing; 0 if fractional
% heuristic: in the case of integral source selection & routing, 1
% indicating that heuristic is used to solve the MMUFP problem ('LP relaxation with randomized rounding'); if 0, solve ILP
% Output:
% cost: total routing cost of the solution
% load_max: load factor (load/capacity) for the most congested link
% flows: P*n matrix; flows(p,i): absolute rate of commodity i on path
% paths{p,i} (P: maximum path index per commodity)
% paths: P*n cell array; paths{p,i}: the node sequence on the p-th path for
% commodity i
% request_type: n*2 matrix; request_type(i,:) denotes that commodity i is
% for serving content item request_type(i,1) to node request_type(i,2). 

% generate auxiliary topology:
[C,V] = size(request); % #items and #nodes
G = [G zeros(V,C); zeros(C,V+C)];
c_link = [c_link zeros(V,C); zeros(C,V+C)];
total_request = sum(sum(request));
for i=1:C
    G(V+i,placement(:,i)==1) = 1; % virtual link (v_i, v) for each v storing i, with cost 1 (to be removed later)
    c_link(V+i,placement(:,i)==1) = total_request; % unlimited capacity on virtual links
end
n = nnz(request); % #commodities (types of requests with nonzero rate)
request_type = zeros(n,2); % request_type(i,:) describes type-i request
s = zeros(n,1);
d = zeros(n,1);
lambda = zeros(n,1);
j = 0;
for i=1:C
    for v=1:V
        if request(i,v) > 0
            j = j + 1;
            request_type(j,:) = [i v];
            lambda(j) = request(i,v);
            s(j) = V+i;
            d(j) = v;
        end
    end
end
E = nnz(G); % #directed links (including virtual links)
links = zeros(E,2); % links(j,:): the two endpoints of link j
link_index = zeros(size(G)); % link_index(u,v): index of link (u,v)
j = 0;
for u=1:length(G)
    for v=1:length(G)
        if G(u,v) ~= 0 
            j = j + 1;
            links(j,:) = [u v];
            link_index(u,v) = j;
        end
    end
end

% solving minimum-cost multiple-source splittable/unsplittable flow
% problem:
n_large = 10^6; % threshold for large matrix/vector
V1 = V+C; % #nodes in the extended G
if max(E^2*n, V1*E*n^2) <= n_large
    A = zeros(E,E*n);
    Aeq = zeros(V1*n,E*n);
    beq = zeros(V1*n,1);
else
    A = sparse(E,E*n);
    Aeq = sparse(V1*n,E*n);
    beq = sparse(V1*n,1);
end
f = zeros(E*n,1);
for i=1:n
    for e=1:E
        f((i-1)*E+e) = lambda(i)*G(links(e,1),links(e,2));
    end
end% min f'*X
b = zeros(E,1);
for e=1:E
    A(e,([1:n]-1)*E+e) = reshape(lambda,1,n);
    b(e) = c_link(links(e,1),links(e,2));
end % A*X <= b (link capacity constraint)
for i=1:n
    for u=1:V1
        for w=1:V1
            if link_index(u,w) > 0
                Aeq((i-1)*V1+u, (i-1)*E+link_index(u,w)) = 1;
            end
            if link_index(w,u) > 0
                Aeq((i-1)*V1+u, (i-1)*E+link_index(w,u)) = -1;
            end
        end
    end
end
for i=1:n
    beq((i-1)*V1+s(i)) = 1;
    beq((i-1)*V1+d(i)) = -1;
end% Aeq*X = beq (flow conservation constraint)
lb = zeros(E*n,1);
ub = ones(E*n,1);

if integral && ~heuristic % ILP
    intcon = 1:E*n;
    options = optimoptions('intlinprog','Display','none');
    [X,cost,~,~] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub, options);
    var_f = reshape(X,E,n); % integral link-level flow
else% ~integral || heuristic % LP
    options = optimoptions('linprog','Display','none');
    [X,cost] = linprog(f,A,b,Aeq,beq,lb,ub,options);
    var_f = reshape(X,E,n); % var_f(e,i): fraction of commodity i routed on link e (fractional link-level flow)
end
var_f_path = zeros(E,n); % var_f_path(p,i): fraction of commodity i routed on p-th path
paths = cell(E,n); % paths{p,i}: p-th path for commodity i as a node sequence
costs = zeros(E,n); % costs(p,i): sum of link costs on p-th path for commodity i
for i=1:n
    A = zeros(V1);
    for e=1:E
        A(links(e,1),links(e,2)) = var_f(e,i);
    end
    [dist, dt, ft, pre] = dfs(sparse(A),s(i));
    p = 0;
    while dt(d(i)) ~= -1 % find a new path
        p = p + 1;
        paths{p,i} = [pre(d(i)) d(i)];
        costs(p,i) = G(pre(d(i)), d(i));
        var_f_path(p,i) = A(pre(d(i)),d(i));
        v = d(i);
        while pre(v) ~= s(i)
            v = pre(v);
            paths{p,i} = [pre(v) paths{p,i}];
            costs(p,i) = costs(p,i) + G(pre(v),v);
            var_f_path(p,i) = min(var_f_path(p,i),A(pre(v),v));
        end
        for j=1:length(paths{p,i})-1
            A(paths{p,i}(j),paths{p,i}(j+1)) = A(paths{p,i}(j),paths{p,i}(j+1)) - var_f_path(p,i);
        end
        [dist, dt, ft, pre] = dfs(sparse(A),s(i));
    end
end
if integral && heuristic % LP relaxation: need to round to integral solution
    cost = 0;
    for i=1:n
        p = find(rand<=cumsum(var_f_path(:,i)),1,'first'); % select each path p with probability var_f_path(p,i)
        var_f_path(:,i) = 0;
        var_f_path(p,i) = 1;
        cost = cost + lambda(i)*costs(p,i);
    end
end% path-level flow recorded in 'var_f_path' and 'paths', with total routing cost 'cost'
p_max = 0; % maximum #paths per commodity
for i=1:n
    p_max = max(p_max,find(var_f_path(:,i)>0,1,'last'));
end
var_f_path = var_f_path(1:p_max,:);
paths = paths(1:p_max,:); % ignore the unused paths (to save space and time)

% convert solution back to original topology:
cost = cost - total_request; % ignore the cost on virtual links
flows = zeros(size(var_f_path)); % flows(p,i): absolute rate of commodity i on path p
load = zeros(V); % load(u,v): load factor for link (u,v)
for i=1:n
    for p=1:p_max
        if ~isempty(paths{p,i})
            paths{p,i} = paths{p,i}(2:end); % remove the virtual source from each path
        end
        flows(p,i) = var_f_path(p,i) * lambda(i); % absolute rate of commodity i on path p
        for k=1:length(paths{p,i})-1
            load(paths{p,i}(k),paths{p,i}(k+1)) = load(paths{p,i}(k),paths{p,i}(k+1)) + flows(p,i)/c_link(paths{p,i}(k),paths{p,i}(k+1)); 
        end
    end
end
load_max = max(max(load)); 

end

