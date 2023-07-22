function [c_RNR, var_x, var_r ] = integral_caching_RNR_Ioannidis( G, c_cache, request, candidatePaths, client )

%Integral caching and integral routing under unlimited link capacities
%according to "offline source routing" in [Ioannidis18JSAC]
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_cache: V*1 vector of cache capacities; c_cache(v): capacity of cache at
% node v. 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% candidatePaths: candidatePaths{i}{j}: j-th shortest path from node i (a client) to the server
% client: node indices of the client nodes (requesters)
% Output:
% c_RNR: routing cost under the given solution
% var_x: V*C matrix, var_x(v,i) = x_{vi}
% var_r: K*Vc*C matrix, var_r(p,s,i) = r_{(i,client(s)),p}
[C,V] = size(request); % #items and #nodes
Vc = length(client); % #clients
K = length(candidatePaths{client(1)}); % #candidate paths per client
for s=2:length(client)
    K = max(K, length(candidatePaths{client(s)}));
end
pmax = 0; % max(path length in #nodes) - 1
for i=1:length(client)
    for p=1:length(candidatePaths{client(i)})
        pmax = max(length(candidatePaths{client(i)}{p})-1, pmax);
    end
end
%% solve LP for L_SR maximization:
% X = [x; r; z]
n_large = 10^6; % threshold for large matrix/vector
if (V*C + K*Vc*C + C*Vc*K*pmax)*max(V,C*Vc*K*pmax) <= n_large
    f = zeros(V*C + K*Vc*C + C*Vc*K*pmax,1);
    A1 = zeros(C*Vc*K*pmax, length(f));
    A2 = zeros(V, length(f));
    A3 = zeros(C*Vc, length(f));
else
    f = sparse(V*C + K*Vc*C + C*Vc*K*pmax,1);
    A1 = sparse(C*Vc*K*pmax, length(f));
    A2 = sparse(V, length(f));
    A3 = sparse(C*Vc, length(f));
end
for i=1:C
    for s=1:Vc
        for p=1:length(candidatePaths{client(s)})% K
            for k=1:length(candidatePaths{client(s)}{p})-1% pmax
                f(V*C + K*Vc*C + (i-1)*Vc*K*pmax+(s-1)*K*pmax+(p-1)*pmax+k ) = - request(i,client(s)) * G(candidatePaths{client(s)}{p}(k), candidatePaths{client(s)}{p}(k+1));
            end
        end
    end
end% obj: min f'*X
for i=1:C
    for s=1:Vc
        for p=1:length(candidatePaths{client(s)})% K
            for k=1:length(candidatePaths{client(s)}{p})-1% pmax
                A1((i-1)*Vc*K*pmax+(s-1)*K*pmax+(p-1)*pmax+k, (i-1)*V+ candidatePaths{client(s)}{p}(k+1: end) ) = -1;  
                A1((i-1)*Vc*K*pmax+(s-1)*K*pmax+(p-1)*pmax+k, V*C + (i-1)*Vc*K+(s-1)*K+p) = 1;
                A1((i-1)*Vc*K*pmax+(s-1)*K*pmax+(p-1)*pmax+k, V*C + K*Vc*C + (i-1)*Vc*K*pmax+(s-1)*K*pmax+(p-1)*pmax+k ) = 1;
            end
        end
    end
end
b1 = ones(C*Vc*K*pmax,1); % A1*X<=b1
for v=1:V
    A2(v, ([1:C]-1)*V+v) = 1;
end
b2 = reshape(c_cache,V,1); % A2*X = b2
for i=1:C
    for s=1:Vc
        A3((i-1)*Vc+s, V*C + (i-1)*Vc*K+(s-1)*K+[1:length(candidatePaths{client(s)})]) = 1;
    end
end
b3 = ones(C*Vc,1); % A3*X = b3
lb = [zeros(V*C+K*Vc*C,1); -inf(C*Vc*K*pmax,1)];
ub = ones(size(f));
options = optimoptions('linprog','Display','none');
[X, fval] = linprog(f, A1,b1, [A2; A3],[b2; b3], lb, ub, options);
var_x = reshape(X(1:V*C),V,C); % var_x(v,i) = x_{vi}
var_r = reshape(X(V*C + [1:K*Vc*C]), K, Vc, C); % var_r(p,s,i) = r_{(i,client(s)),p}

frac_x = var_x; 
frac_r = var_r;
%% round var_x to integral caching:
I = [];
for v=1:V
    if any(var_x(v,:)>0 & var_x(v,:)<1)
        I = find(var_x(v,:)>0 & var_x(v,:)<1,2); % found two fractional placements var_x(v,I(1)), var_x(v,I(2))
        break;
    end
end
while ~isempty(I) && length(I)==2
    x_candidate = cell(1,2); 
    for i=1:2
        j = setdiff(1:2,i);
        x_candidate{i} = var_x;
        delta = min(var_x(v,I(i)), 1-var_x(v,I(j)));
        x_candidate{i}(v,I(i)) = x_candidate{i}(v,I(i)) - delta;
        x_candidate{i}(v,I(j)) = x_candidate{i}(v,I(j)) + delta;
    end
    if C_RNR_Ioannidis( x_candidate{1}, var_r, request, candidatePaths, G, client ) < C_RNR_Ioannidis( x_candidate{2}, var_r, request, candidatePaths, G, client )
        var_x = x_candidate{1};
    else
        var_x = x_candidate{2};
    end
    I=[];
    for v=1:V
        if any(var_x(v,:)>0 & var_x(v,:)<1)
            I = find(var_x(v,:)>0 & var_x(v,:)<1,2);
            break;
        end
    end
end% round var_x to an integral solution
%% update var_r to RNR:
c_RNR = 0;
for i=1:C
    for s=1:length(client)
        pmin = 0; cmin = inf; 
        for p = 1:length(candidatePaths{client(s)})
            c = 0;
            for k=1:length(candidatePaths{client(s)}{p})-1
                c = c + G(candidatePaths{client(s)}{p}(k),candidatePaths{client(s)}{p}(k+1))*prod(1-var_x(candidatePaths{client(s)}{p}(k+1:end),i));
            end
            if c<cmin
                pmin = p; cmin = c;
            end
        end
        c_RNR = c_RNR + request(i,client(s))*cmin; 
        var_r(:,s,i) = 0; var_r(pmin,s,i) = 1;
    end
end

%% sanity check:
if false
    
% LP solver correctness:
l = L_RNR_Ioannidis( var_x, var_r, request, candidatePaths, G, client );
l_frac = L_RNR_Ioannidis( frac_x, frac_r, request, candidatePaths, G, client );
if abs(l_frac + fval) > 10^(-6)
    error(['LP objective should equal -L_SR(frac_x, frac_r)']);
end
if l_frac < l
    error('l_frac should be no smaller than l');
end
% rounding of var_x correctness:
c1 = C_RNR_Ioannidis( var_x, frac_r, request, candidatePaths, G, client );
c2 = C_RNR_Ioannidis( frac_x, frac_r, request, candidatePaths, G, client );
if c1 > c2
    error('rounding of var_x should not increase routing cost');
end
% correctness of Lemma 2 in [Ioannidis18JSAC]:
c0 = 0;
for i=1:C
    for s=1:length(client)
        for p=1:length(candidatePaths{client(s)})
            for k=1:length(candidatePaths{client(s)}{p})-1
                c0 = c0 + request(i,client(s))*G(candidatePaths{client(s)}{p}(k),candidatePaths{client(s)}{p}(k+1));
            end
        end
    end
end
f = c0-c2;
if f < (1-1/exp(1))*l_frac-10^(-6) || f > l_frac + 10^(-6)
    error(['Lemma 2 in [Ioannidis18JSAC] is violated: f = ' num2str(f) ', l_frac = ' num2str(l_frac)]);
end
  
end



end

