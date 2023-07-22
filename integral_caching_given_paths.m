function [ cost, load_max, var_x ] = integral_caching_given_paths( G, c_cache, c_link, request, flows, paths, request_type )
%(1-1/e)-approximation for integral content placement under given path(s)
%for each commodity (possibly split over multiple paths)
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_cache: V*1 vector of cache capacities; c_cache(v): capacity of cache at
% node v. 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% flows: P*n matrix; flows(p,i): absolute rate of commodity i on path p
% paths: P*n cell array; paths{p,i}: the node sequence on path p for
% commodity i
% request_type: n*2 matrix; request_type(i,:) denotes that commodity i is
% for serving content item request_type(i,1) to node request_type(i,2). 
% Output:
% cost: total routing cost under the solution
% var_x: V*C matrix; var_x(v,i) = x_{vi} (whether node v stores item i)

[C,V] = size(request); % #items and #nodes
[P,n] = size(flows); % #paths per commodity, #commodities
%% solve LP for L_{r,f} maximization
% X = [x; z]
n_large = 10^6; % threshold for large matrix/vector
if n*P*V*(V*C+n*P*V) <= n_large
    A = zeros(n*P*V,V*C+n*P*V);
    Aeq = zeros(V, V*C+n*P*V);
else
    A = sparse(n*P*V,V*C+n*P*V);
    Aeq = sparse(V, V*C+n*P*V);
end
f = zeros(V*C+n*P*V,1);
for i=1:n
    for p=1:P
        for k=1:length(paths{p,i})-1
            path = paths{p,i}; 
            f(V*C+(i-1)*P*V+(p-1)*V+k) = - flows(p,i)*G(path(length(path)-k),path(length(path)-k+1));
        end
    end
end% min f'*X
for i=1:n
    for p=1:P
        for k=1:length(paths{p,i})-1
            A((i-1)*P*V+(p-1)*V+k, (request_type(i,1)-1)*V+paths{p,i}(length(paths{p,i})- [0:k-1]) ) = -1;
            A((i-1)*P*V+(p-1)*V+k, V*C+ (i-1)*P*V+(p-1)*V+k) = 1;
        end
    end
end
b = zeros(n*P*V,1); % part of the definition: z^i_{p,k} <= sum_{k'=0}^{k-1}x_{p(|p|-k'),request_type(i,1)}
for v=1:V
    Aeq(v, ([1:C]-1)*V+v) = 1;
end
beq = reshape(c_cache,V,1); % Aeq = beq (cache capacity constraint; assumed to use up capacity to enable rounding)
ub = ones(V*C+n*P*V,1);
lb = [zeros(V*C,1); -inf(n*P*V,1)];
for i=1:n
    for p=1:P
        if flows(p,i)>0
            lb((request_type(i,1)-1)*V + paths{p,i}(1)) = 1; %node 'paths{p,i}(1)' must have content 'request_type(i,1)'
        end
    end
end% the source of a path used to serve (possibly part of) a content must store that content. 
options = optimoptions('linprog','Display','none');
[X, fval] = linprog(f,A,b,Aeq,beq,lb,ub,options);
var_x = reshape(X(1:V*C),V,C); % var_x(v,i) = x_{vi}

frac_x = var_x; 
%% round fractional solution to integral solution:
I=[];
for v=1:V
    if any(var_x(v,:)>0 & var_x(v,:)<1)
        I = find(var_x(v,:)>0 & var_x(v,:)<1,2);
        break;
    end
end
while ~isempty(I) && length(I)==2 % var_x(v,I(1)) and var_x(v,I(2)) \in (0,1):  
    coeff = zeros(1,2); % coefficient of var_x(v,I) in F_{r,f}
    for j=1:2 % for var_x(v,I(j))
        commodity = find(request_type(:,1)==I(j));
        for i = commodity % commodities requesting item I(j)
            for p=1:P
                if any(paths{p,i}==v)
                    kv = length(paths{p,i}) - find(paths{p,i}==v);
                    for k=kv+1:length(paths{p,i})-1
                        coeff(j) = coeff(j) + flows(p,i)*G(paths{p,i}(length(paths{p,i})-k),paths{p,i}(length(paths{p,i})-k+1))*prod(1-var_x(paths{p,i}(length(paths{p,i})-[0:kv-1 kv+1:k-1]) ,I(j)));
                    end
                end
            end
        end
    end
    [~,i_max]=max(coeff); i_min = setdiff(1:2,i_max);
    x_vi = min(1,var_x(v,I(1))+var_x(v,I(2))); 
    x_vj = var_x(v,I(1))+var_x(v,I(2)) - x_vi;
    var_x(v,I(i_max)) = x_vi;
    var_x(v,I(i_min)) = x_vj;
    I=[];
    for v=1:V
        if any(var_x(v,:)>0 & var_x(v,:)<1)
            I = find(var_x(v,:)>0 & var_x(v,:)<1,2);
            break;
        end
    end
end% round var_x to an integral solution
%% compute total routing cost under this solution:
cost = C_r_f(G, flows, paths, request_type, var_x);

load = zeros(V); % load(u,v): load factor for link (u,v)
for i=1:n
    for p=1:P
        if ~isempty(paths{p,i})
%             disp(['integral_caching_given_paths: trim paths{' num2str(p) ',' num2str(i) '}'])
            s = find(var_x(paths{p,i},request_type(i,1))>0,1,'last'); % paths{p,i}(s) is the node closest to destination on paths{p,i} that stores content i
            paths{p,i} = paths{p,i}(s:end); % trim the path to start from the cache storing content i that is closest to the requester
            for k=1:length(paths{p,i})-1
                load(paths{p,i}(k),paths{p,i}(k+1)) = load(paths{p,i}(k),paths{p,i}(k+1)) + flows(p,i)/c_link(paths{p,i}(k),paths{p,i}(k+1)); 
            end
        end
    end
end
load_max = max(max(load)); 

if false
%% sanity check:
L_int = L_r_f(G, flows, paths, request_type, var_x);
L_frac = L_r_f(G, flows, paths, request_type, frac_x);
if abs(L_frac + fval) > 10^(-10) || L_int > L_frac
    error(['frac_x did not maximize L_{r,f}(x)']);
end
cost_frac = C_r_f(G, flows, paths, request_type, frac_x);
if cost_frac < cost
    error(['rounded cost ' num2str(cost) ' should not be larger than fractional cost ' num2str(cost_frac)]);
end
% c0 = 0;
% for i=1:n
%     for p=1:P
%         for k=1:length(paths{p,i})-1
%             c0 = c0 + flows(p,i)*G(paths{p,i}(length(paths{p,i})-k),paths{p,i}(length(paths{p,i})-k+1));
%         end
%     end
% end% F_r_f = c0 - C_r_f
end

end