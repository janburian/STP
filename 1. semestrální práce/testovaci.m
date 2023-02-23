%% https://github.com/kratochviljakub

%% 1
clear
P = [0.2 0.2 0.2 0.2 0.1 0.1;
     0.3 0.1 0.1 0.1 0.1 0.3;
     0.5 0.1 0.1 0.1 0.1 0.1;
     0.1 0.2 0.2 0.2 0.2 0.1;
     0.2 0.1 0.1 0.3 0.1 0.2;
     0.1 0.2 0.3 0.1 0.1 0.2];

mc = dtmc(P);
graphplot(mc,'ColorEdges',true)
n = rank(P);
I = eye(n);
B = [[I-P]'; ones(1,n)];
f = [zeros(n,1); 1];
a = (B\f)'
A = [a; a; a; a; a; a];
Z = inv(I - (P - A));
Z1 = diag(diag(Z));
M1 = inv(diag(diag(A)));
E = ones(n);
M = (I - Z + E * Z1) * M1

 %% 2
clear
P = [1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
     0.0 0.4 0.0 0.0 0.4 0.2;
     0.6 0.0 0.0 0.2 0.2 0.0;
     0.2 0.2 0.2 0.2 0.0 0.2;
     0.1 0.1 0.1 0.1 0.1 0.5];
 
mc = dtmc(P);
graphplot(mc,'ColorNodes',true,'ColorEdges',true)
Q = P(3:6,3:6);
T = inv(eye(4) - Q)
t = T * ones(4,1)
R = P(3:6,1:2);
d = T * R

%%
A = [ 0.1 -0.3  0.3; 0.1 0.2 -0.35; 1 1  1];
B = [0; 0; 1];
X = linsolve(A,B)