function [ cost ] = L_r_f( G, flows, paths, request_type, var_x )
%Concave relaxation L_{r,f}(x):
cost = 0;
[P,n] = size(flows);
for i=1:n
    for p=1:P
        for k=1:length(paths{p,i})-1
            cost = cost + flows(p,i)*G(paths{p,i}(length(paths{p,i})-k), paths{p,i}(length(paths{p,i})-k+1))*min(1,sum(var_x(paths{p,i}(length(paths{p,i})-[0:k-1]), request_type(i,1))));
        end
    end
end


end

