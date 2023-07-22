% toy example for debugging:
G = zeros(5);
G(2,1) = 1;
G(3,1)=2;
G(4,1)=3;
G(5,2)=1;
G(5,3)=1;
G(5,4)=1; % set virtual link costs = 1 (as G(.,.)=0 is interpreted as no link)
% the final cost needs to be reduced by sum(sum(request))*1 to ignore this
% cost of virtual links. 
s=5;
d=[1 1 1];
lambda = [.9,.8,.7];
c_link = zeros(size(G));
c_link(2,1)=1;
c_link(3,1)=1;
c_link(4,1)=1;
c_link(5,2)=sum(lambda); % effectively infinite capacity
c_link(5,3)=c_link(5,2);
c_link(5,4)=c_link(5,2);

Kmax = 10;
Cost = zeros(1,Kmax);
Load = zeros(1,Kmax);
for K=1:Kmax
    disp(['K = ' num2str(K) ':'])
    [ Cost(K), Load(K), paths_output ] = bicriteria_MSUFP( G, c_link, s, d, lambda, K );
end
figure;
subplot(2,1,1)
plot(1:Kmax,Cost);
xlabel('K')
ylabel('cost')
subplot(2,1,2)
plot(1:Kmax,Load);
xlabel('K')
ylabel('load factor')
