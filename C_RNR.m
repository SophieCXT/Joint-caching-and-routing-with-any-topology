function [ val ] = C_RNR( var_x, var_r, C, V, request, SPcost, w_max )
%C_RNR(x,r) as defined in (2):
% total routing cost under content placement x, source selection r, and
% shortest path routing.
val = 0;
for i=1:C
    for s=1:V
        val = val + request(i,s)*sum(var_r(:,s,i).*(var_x(:,i).*SPcost(:,s) + (1-var_x(:,i))*w_max));
    end
end

end

