function test_gurobi_LP()
addpath("C:\gurobi912\win64\matlab");
% >> cd C:\gurobi912\win64\matlab
% >> gurobi_setup
%% ?????????
% maxf(x)=2*x1+3*x2
%    1       x1+2*x2<=3
%    2       4*x1   <=16
%    3          4*x2<=12
%    4       x1      >=0
%    5             x2>=0
%     A=[1 2;4 0;0 4;1 0;0 1];
%     obj=[2 3]
%     rhs=[3 16 12 0 0]
%     sense=['<' '<' '<' '>' '>']
model.A=sparse([1 2;4 0;0 4;1 0;0 1]);
model.obj=[2 3];
model.modelsense='MAX';
model.rhs=[3 16 12 0 0];
model.sense=['<' '<' '<' '>' '>'];
 
% ????????
result = gurobi(model);
disp(result.objval);
disp(result.x);
% ????
model.vbasis = result.vbasis;
model.cbasis = result.cbasis;
result = gurobi(model);
end

