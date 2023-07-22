% toy example for debugging:
G = zeros(5);
G(2,4)=rand;
G(4,1)=rand;
G(3,4)=50;
G(4,5)=100;
G=max(G,G');
c_cache=[0 1 1 0 2]';
request=zeros(2,5);
request(1,1)=2;
request(2,1)=1;

[ c_RNR, var_x, var_r ] = integral_caching_RNR( G, c_cache, request )
c_sub = request(1,1)*(G(3,4)+G(4,1))+request(2,1)*(G(2,4)+G(4,1))
c_opt = request(1,1)*(G(2,4)+G(4,1))+request(2,1)*(G(3,4)+G(4,1))

var_x_opt = [0 0; 1 0; 0 1; 0 0; 1 1];
for i=1:C
    for s=1:V
        var_r_opt(:,s,i) = 0;
        [c_is,v_is] = min(var_x_opt(:,i).*SPcost(:,s) + (1-var_x_opt(:,i))*w_max);
        var_r_opt(v_is,s,i) = 1;% v_is has the nearest replica for request (i,s)
    end
end
l_RNR_opt = L_RNR( var_x_opt, var_r_opt, C, V, request, w_max, SPcost )