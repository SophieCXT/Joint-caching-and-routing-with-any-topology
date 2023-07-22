function [ c_RNR, var_x, var_r ] = integral_caching_RNR( G, c_cache, request )
%Integral caching and integral routing under unlimited link capacities
%(Algorithm 1)
% Input:
% G: V*V adjacency matrix containing link costs; G(i,j) = 0 means link (i,j)
% does not exist; G(i,j) > 0 means link (i,j) has cost G(i,j) (assuming
% positive costs on all existing links). 
% c_cache: V*1 vector of cache capacities; c_cache(v): capacity of cache at
% node v. 
% request: C*V demand matrix; request(i,s): request rate for item i
% by node s.
% Output:
% c_RNR: routing cost under the given solution
% var_x: V*C matrix for integral caching decision (var_x(v,i) = x_{vi}: whether node v
% stores item i)
% var_r: integral source selection decision (var_r(v,s,i) = r^{(i,s)}_v: whether request (i,s) is served by node v)

[C,V] = size(request); % #items and #nodes
SP = cell(1,V); % SP{s}{d}: least-cost path as a node sequence from s to d
SPcost = zeros(V); % SPcost(s,d): least cost from s to d

%% line 1:
for s=1:V
    [SPcost(s,:), SP{s}] = Dijkstra_source(G, s); 
end
w_max = max(max(SPcost)) + 1; 

%% line 2:
% X = [x; r; z]: column-vector of decision variables
n_large = 10^6; % threshold for large matrix/vector
if C*V^2*(C*V+2*C*V^2) <= n_large
    f=zeros(V*C+2*V^2*C,1);
    A1 = [zeros(C*V^2,C*V+C*V^2) eye(C*V^2)];
    A2 = [zeros(C*V^2,C*V) eye(C*V^2) eye(C*V^2)];
    A3 = zeros(V, C*V+2*C*V^2);
    A4 = zeros(C*V, C*V+2*C*V^2);
else
    f=sparse(V*C+2*V^2*C,1);
    A_eye = sparse(C*V^2,C*V^2);
    for i=1:C*V^2
        A_eye(i,i) = 1;
    end% A_eye = sparse(eye(C*V^2))    
    A1 = [sparse(C*V^2,C*V+C*V^2) A_eye];
    A2 = [sparse(C*V^2,C*V) A_eye A_eye];
    A3 = sparse(V, C*V+2*C*V^2);
    A4 = sparse(C*V, C*V+2*C*V^2);
end
for i=1:C
    for s=1:V
        f(V*C + V^2*C + (i-1)*V^2+(s-1)*V+ [1:V],1) = - request(i,s)*w_max; 
    end 
end% obj: min f'*X (8a)
b1 = ones(C*V^2,1); % A1*X <= b1 (8b)
for i=1:C
    for v=1:V
        A2((i-1)*V^2 + ([1:V]-1)*V + v, (i-1)*V+v) = ((SPcost(v,:)-w_max)./w_max)';
    end
end
b2 = b1; % A2*X <= b2 (8c)
for v=1:V
    A3(v,([1:C]-1)*V+v) = 1;
end
b3 = reshape(c_cache,V,1); % A3*X <= b3 (3c) --> A3*X = b3 (use up cache capacity at every node, to enable rounding)
for i=1:C
    for s=1:V
        A4((i-1)*V+s, C*V + (i-1)*V^2+(s-1)*V+[1:V]) = 1;
    end
end
b4 = ones(C*V,1); % A4*X = b4 (3b)
lb = [zeros(C*V+C*V^2,1); -inf(C*V^2,1)];
ub = [ones(C*V+C*V^2,1); inf(C*V^2,1)]; % lb <= X <= ub (8e)
options = optimoptions('linprog','Display','none');
[X, fval] = linprog(f,[A1; A2],[b1; b2],[A3; A4],[b3; b4],lb,ub,options);
var_x = reshape(X(1:C*V),V,C); % var_x(v,i) = x_{vi}
var_r = reshape(X(C*V+[1:C*V^2]), V,V,C); % var_r(v,s,i) = r^{(i,s)}_v

frac_x = var_x;
frac_r = var_r;

%% line 3:
I=[];
for v=1:V
    if any(var_x(v,:)>0 & var_x(v,:)<1)
        I = find(var_x(v,:)>0 & var_x(v,:)<1,2);
        break;
    end
end
while ~isempty(I) && length(I)==2 % var_x(v,I(1)) and var_x(v,I(2)) \in (0,1):    
    coeff_i = sum(request(I(1),1:V).*var_r(v,1:V,I(1)).*(w_max-SPcost(v,1:V)));
    coeff_j = sum(request(I(2),1:V).*var_r(v,1:V,I(2)).*(w_max-SPcost(v,1:V)));
    if coeff_i >= coeff_j
        x_vi = min(1,var_x(v,I(1))+var_x(v,I(2))); 
        x_vj = var_x(v,I(1))+var_x(v,I(2)) - x_vi;
    else
        x_vj = min(1,var_x(v,I(1))+var_x(v,I(2)));
        x_vi = var_x(v,I(1))+var_x(v,I(2)) - x_vj;
    end
    var_x(v,I(1)) = x_vi;
    var_x(v,I(2)) = x_vj;
    I=[];
    for v=1:V
        if any(var_x(v,:)>0 & var_x(v,:)<1)
            I = find(var_x(v,:)>0 & var_x(v,:)<1,2);
            break;
        end
    end
end% round var_x to an integral solution
%% line 4:
c_RNR = 0; 
for i=1:C
    for s=1:V
        var_r(:,s,i) = 0;
        [c_is,v_is] = min(var_x(:,i).*SPcost(:,s) + (1-var_x(:,i))*w_max);
        var_r(v_is,s,i) = 1;% v_is has the nearest replica for request (i,s)
        c_RNR = c_RNR + request(i,s)*c_is; 
    end
end% var_r is an integral solution by RNR



%% sanity check: (can be commented out later for speed)
if false
% LP solver correctness:
l_RNR = L_RNR( var_x, var_r, C, V, request, w_max, SPcost );
l_RNR_frac = L_RNR( frac_x, frac_r, C, V, request, w_max, SPcost );
if abs(l_RNR_frac + fval) > 10^(-6)
    error(['LP objective should equal -L_RNR(frac_x, frac_r)']);
end
if l_RNR_frac < l_RNR
    error(['L_RNR_frac is no smaller than L_RNR']);
end
% rounding of var_x correctness:
c_RNR1 = C_RNR( var_x, frac_r, C, V, request, SPcost, w_max );
c_RNR2 = C_RNR( frac_x, frac_r, C, V, request, SPcost, w_max );
if c_RNR1 > c_RNR2
    error(['C_RNR(integral_x, fractional_r) is no greater than C_RNR(fractional_x, fractional_r)']);
end
% correctness of Lemma IV.1:
c0_RNR = V*w_max*sum(sum(request));
f_RNR_frac = c0_RNR-c_RNR2; 
if f_RNR_frac < (1-1/exp(1))*l_RNR_frac || f_RNR_frac > l_RNR_frac
    error(['L_RNR and F_RNR need to satisfy the bounds in Lemma IV.1']);
end
end

end

