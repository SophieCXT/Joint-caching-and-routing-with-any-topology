function [ paths_output ] = rounding_MSUFP( G, s, d, lambda_bar, var_f_link, links, link_index  )
%implementation of Algorithm 2 in [Skutella2002]:
% Input:
% G: link costs (G(i,j)=0 means no link from node i to node j)
% s: source node
% d: d(i) is destination of commodity i
% lambda: lambda(i) is demand of commodity i; assume that for any i,j,
% lambda(i) and lambda(j) differ by a factor of 2^q for some integer q
% var_f_link: a splittable flow f_0 satisfying all demands; var_f_link(u,v)
% is the total flow on link (u,v). 
% Output:
% paths_output: cell array, paths_output{i} is the path for commodity i as
% node sequence. 
V = length(G);
n = length(d);
E = nnz(G);
paths_output = cell(1,n);
if n == 0
    return;
end
[~,I] = sort(lambda_bar,'ascend');
i = 1;
j = 0;
A = var_f_link; 
    f = zeros(E,1);
    for e=1:E
        f(e) = G(links(e,1),links(e,2));
    end
    Aeq = zeros(V,E);
    for u=1:V
        for w=1:V
            if link_index(u,w) > 0
                Aeq(u,link_index(u,w)) = 1;
            end
            if link_index(w,u) > 0
                Aeq(u,link_index(w,u)) = -1;
            end
        end
    end
%%
while lambda_bar(I(1))*2^j <= lambda_bar(I(n)) + 10^(-8)
    j = j + 1;
    delta_j = lambda_bar(I(1))*2^(j-1);
    % disp(['j = ' num2str(j) ', delta_j = ' num2str(delta_j)])
    c_link = delta_j.*ceil(A./delta_j); % c_link(u,v): "capacity" of link (u,v)
    % compute a delta_j-integral flow satisfying demands 1:n with cost no
    % worse than var_f_link
    % Method 1: LP (issue: may not find a delta_j-integral solution, because linprog does not use simplex method) 
%     f = zeros(E,1);
%     for e=1:E
%         f(e) = G(links(e,1),links(e,2));
%     end
%     Aeq = zeros(V,E);
%     for u=1:V
%         for w=1:V
%             if link_index(u,w) > 0
%                 Aeq(u,link_index(u,w)) = 1;
%             end
%             if link_index(w,u) > 0
%                 Aeq(u,link_index(w,u)) = -1;
%             end
%         end
%     end
    beq = zeros(V,1);
    beq(s) = sum(lambda_bar(I(i:end))); %sum(lambda_bar);
    for i1=i:length(I)
        beq(d(I(i1))) = beq(d(I(i1))) - lambda_bar(I(i1));
    end% Aeq*x = beq (flow conservation constraint) -> Aeq*x <= beq + epsilon, -Aeq*x<=-beq+epsilon; note: ignore commodities I(1),...I(i-1) that have been routed
%     for u=[1:s-1 s+1:V]
%         beq(u) = -sum(lambda_bar(d==u));
%     end
    lb = zeros(E,1); 
    ub = zeros(E,1);
    for e=1:E
        ub(e) = c_link(links(e,1),links(e,2));
    end
%     epsilon = ones(V,1)*10^(-8);
    x = zeros(E,1);
    for e=1:E
        x(e) = A(links(e,1),links(e,2));
    end% x is the given solution according to the input fractional flow
    epsilon = ones(V,1)*max(abs(Aeq*x-beq));    
    % disp(['rounding_MSUFP: given error margin ' num2str(epsilon(1)) ' in the LP for computing the min-cost delta_j-integral flow'])
%     [f_j,fval] = linprog(f,[],[],Aeq,beq,lb,ub); % var_f_link is a feasible solution to this LP, hence the optimal solution f_j has a cost that is no larger
    options = optimoptions('linprog','Display','none');
    [f_j,fval] = linprog(f,[Aeq; -Aeq],[beq+epsilon; -beq+epsilon],[],[],lb,ub,options); % to tolerate some numerical error
    % Method 2: iterative shortest path routing    
    
    for e=1:E
        if min(f_j(e) - floor(f_j(e)/delta_j)*delta_j, ceil(f_j(e)/delta_j)*delta_j - f_j(e)) > 10^(-6)
            error(['f_j(' num2str(e) ') / delta_j = ' num2str(f_j(e)/delta_j) ', which is not integral']);
        else
            f_j(e) = round(f_j(e)/delta_j)*delta_j; % round to nearest integer multiple of delta_j (to correct numerical error)
        end
    end    
    %A = zeros(V);
    for e=1:E
        A(links(e,1),links(e,2)) = f_j(e);
    end
    while i<=n && abs(lambda_bar(I(i)) - delta_j) < min(10^(-8),delta_j/2) %lambda_bar(I(i)) == delta_j        
        [dist, dt, ft, pre] = dfs(sparse(A>10^(-8)),s);
        paths_output{I(i)} = [d(I(i))];
        v = d(I(i));
        while v ~= s
            paths_output{I(i)} = [pre(v) paths_output{I(i)}];
            % disp(['i = ' num2str(i) ': path_' num2str(I(i)) ' traverses link (' num2str(pre(v)) ', ' num2str(v) '), capacity reduces from ' num2str(A(pre(v),v)) ' to ' num2str(A(pre(v),v) - delta_j)])
            A(pre(v),v) = A(pre(v),v) - delta_j;
            v = pre(v);
        end
        i = i + 1;
    end           
end

if i ~= n+1
    error(['rounding_MSUFP: fail to find paths for all the ' num2str(n) ' commodities']);
end
%% sanity check according to Lemma IV.6:
if false
% (i) cost should be no larger:
cost = 0;
for i=1:n
    for k=1:length(paths_output{i})-1
        cost = cost + lambda_bar(i)*G(paths_output{i}(k),paths_output{i}(k+1));        
    end
end
if cost > sum(sum(G.*var_f_link)) + 10^(-8)
    error(['rounded cost ' num2str(cost) ' is greater than the cost of initial flow ' num2str( sum(sum(G.*var_f_link)) )]);
end
% (ii) link load should be no larger after ignoring the largest flow:
max_flow = zeros(E,1); % max per-commodity flow on each link
total_flow = zeros(E,1); % total flow on each link
for i=1:n
    for k=1:length(paths_output{i})-1
        total_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))) = total_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))) + lambda_bar(i);
        max_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))) = max(max_flow(link_index(paths_output{i}(k),paths_output{i}(k+1))), lambda_bar(i));
    end
end
for e=1:E
%     disp(['link ' num2str(e) ':'])
%     disp(['rounded load: ' num2str(total_flow(e)) ', rounded load without largest demand: ' num2str(total_flow(e) - max_flow(e)) ', initial load: ' num2str(var_f_link(links(e,1),links(e,2)))])
    if var_f_link(links(e,1),links(e,2))>0 && total_flow(e) - max_flow(e) >= var_f_link(links(e,1),links(e,2))
        error(['rounded load (ignoring the largest demand) on link ' num2str(e) ' is ' num2str(total_flow(e) - max_flow(e)) ', not smaller than initial load ' num2str(var_f_link(links(e,1),links(e,2)))]);
    end
end
end

end

