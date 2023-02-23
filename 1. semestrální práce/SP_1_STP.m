%% Prvni semestralni prace z predmetu STP
% Jan Burian

%% Priklad 1
% Markovsky retezec
    % - homogenni (pij(n) nezavisi na n)
    % - regularni (P^(n) pro konecne n (neobsahuje zadne nulove prvky)
% 6 uzlu => matice prechodu P o velikosti 6x6

clc
clear all
close all 

P = [0.45 0.23 0.1 0.05 0.07 0.1; 
     0.35 0.15 0.1 0.2  0.1  0.1;
     0.47 0.03 0.2 0.1  0.05 0.15;
     0.38 0.22 0.3 0.02 0.06 0.02;
     0.23 0.17 0.2 0.12 0.22 0.06;
     0.12 0.08 0.3 0.15 0.3  0.05];
 
mc = dtmc(P);
figure;
graphplot(mc,'ColorEdges',true);

syms a1 a2 a3 a4 a5 a6
a = [a1 a2 a3 a4 a5 a6];
soucin_aP = a*P;
eqns = [soucin_aP(1) == a1,
        soucin_aP(2) == a2,
        soucin_aP(3) == a3,
        soucin_aP(4) == a4,
        soucin_aP(5) == a5,
        soucin_aP(6) == a6,
        a1 + a2 + a3 + a4 + a5 + a6 == 1];
result = solve(eqns, [a1 a2 a3 a4 a5 a6]);
a = structfun(@double,result)'
A = [a; a; a; a; a; a];

n = rank(P); % hodnost matice P 
I = eye(n); % jednotkova matice 
E = ones(n); % matice samych jednicek
Z = inv(I - (P - A));

Z_pomocna = diag(diag(Z));
M_pomocna = inv(diag(diag(A)));

M = (I - Z + E * Z_pomocna) * M_pomocna

%% Priklad 2
% Markovsky retezec
    % - homogenni (pij(n) nezavisi na n)
    % - absorpcni se dvema absorpcnimi stavy
% 6 uzlu => matice prechodu P o velikosti 6x6

clc
clear all
close all

P = [1.0 0.0 0.0 0.0 0.0 0.0; 
     0.0 1.0 0.0 0.0 0.0 0.0;
     0.4 0.3 0.2 0.1 0.0 0.0;
     0.2 0.4 0.2 0.2 0.0 0.0;
     0.1 0.2 0.7 0.0 0.0 0.0;
     0.3 0.1 0.1 0.2 0.2 0.1];

mc = dtmc(P);
graphplot(mc,'ColorNodes',true,'ColorEdges',true)

% kanonicky tvar P = [I 0; R Q]
Q = P(3:6, 3:6); % matice Q vznikne z P, a to vynechanim radku a sloupcu odpovidajicich absorpcnim stavum
I = eye(4); % jednotkova matice o hodnosti 4 = rank(Q)
T = inv(I - Q) % stredni pocet pruchodu stavem j... (fundamentalni matice T)
t = T * ones(4,1) % doba pohybu v tranzientnim stavu
R = P(3:6, 1:2); % matice typu s x r obsahujici pravdepodobnosti prechodu z neabsorpcnich do 
                            % absorpcnich stavu a
                            
d = T * R






