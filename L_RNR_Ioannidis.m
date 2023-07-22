function [ val ] = L_RNR_Ioannidis( var_x, var_r, request, candidatePaths, G, client )
%L_SR(x,r) as defined in (30) of [Ioannidis18JSAC]:
% total routing cost under content placement x, source selection r, and
% shortest path routing.
[C,V] = size(request);
val = 0;
for i=1:C
    for s=1:length(client)
        for p=1:length(candidatePaths{client(s)})
            for k=1:length(candidatePaths{client(s)}{p})-1
                val = val + request(i,client(s))*G(candidatePaths{client(s)}{p}(k), candidatePaths{client(s)}{p}(k+1))*min(1, 1-var_r(p,s,i)+sum(var_x(candidatePaths{client(s)}{p}(k+1:end), i)));
            end
        end
    end
end

end

