function [ val ] = L_RNR( var_x, var_r, C, V, request, w_max, SPcost )
%return L_RNR(x,r) as defined in (5) of Lemma IV.1
%   concave relaxation of cost saving function F_RNR(x,r)
val = 0;
for i=1:C
    for s=1:V
        val = val + request(i,s)*w_max*sum(min(1,1-var_r(:,s,i)+var_x(:,i).*(w_max-SPcost(:,s))./w_max));
    end
end

end

