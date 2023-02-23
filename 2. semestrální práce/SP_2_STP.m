%% Druha semestralni prace z predmetu STP
% Jan Burian

%% Priklad 1 
clc
clear all 
close all 

% Zadane parametry 
Q = 3; 
b = 0.5;
T = 1; 
M = 10^4;
N = 100; 
TAU = 6; % {0, 1, 2, 3, 4, 5}
K = 95; % {0, 1, 2,..., 94}
sigma_bilysum = Q*(1-exp(-2*b*T));

% Generovani realizaci
Xk = zeros(N, M);
for i = 1:M
    Xk(1, i) = randn * sqrt(Q); % vzdy prvni prvek ve sloupci (X0) -> celkem M sloupcu
    for j = 2:N % generovani ostatnich prvku ve sloupci
        % X_{k+1} = e^{-bT}  * X_k          + W_k
        Xk(j, i) = exp(-b*T) * Xk(j - 1, i) + randn * sqrt(sigma_bilysum);
    end
end

% Odhad parametru 
cov_Xk_odhad = zeros(TAU, K);
for k = 1:K
    for tau = 1:TAU
        Xk_x = Xk(k, :); % k-ty casovy okamzik
        Xk_y = Xk(k + tau - 1, :); % k + tau casovy okamzik
        EXk_x = mean(Xk_x); % stredni hodnota
        EXk_y = mean(Xk_y); % stredni hodnota
        for i = 1:length(Xk_x)
            cov_Xk_odhad(tau, k) = cov_Xk_odhad(tau, k) + ((Xk_x(i) - EXk_x) * (Xk_y(i) - EXk_y));
        end
        cov_Xk_odhad(tau, k) = cov_Xk_odhad(tau, k) / length(Xk_x);
    end
end

% Teoreticky vypocet autokovariancni funkce
cov_Xk_vypocet = zeros(TAU, K);
for tau = 1:TAU
    cov_Xk_vypocet(tau, :) = Q * exp((-tau + 1) * b); 
end


% Vykresleni
figure
for j = 1:TAU
    hold on  
    plot(0:length(cov_Xk_odhad)-1, cov_Xk_odhad(j,:), 'LineWidth', 1.5); 
end

hold on
for u = 1:TAU
    hold on
    plot(0:length(cov_Xk_vypocet)-1, cov_Xk_vypocet(u, :), 'black:','LineWidth', 1.5);
end

Legend = cell(7,1);
Legend{1} = 'COV[X_k,X_k]'; 
Legend{2} = 'COV[X_k,X_{k+1}]';  
Legend{3} = 'COV[X_k,X_{k+2}]'; 
Legend{4} = 'COV[X_k,X_{k+3}]';
Legend{5} = 'COV[X_k,X_{k+4}]'; 
Legend{6} = 'COV[X_k,X_{k+5}]'; 
Legend{7} = 'Teoreticke hodnoty \newline jednotlivych COV[X_k,X_{k+\tau} ]'; 

title('Autokovariancni funkce');
xlabel('Cas [k]');
legend(Legend);

% Setreni stacinonarity
figure
for k = 1:K
    Xk_x = Xk(k, :); % k-ty casovy okamzik
    EXk_x = mean(Xk_x); % stredni hodnota
    hold on
    scatter(k, EXk_x)
end
title('Stredni hodnoty casovych okamziku')
xlabel('Cas [k]');

%% Priklad 2
clc
clear all 
close all 

% Zadane parametry 
M = 10^4;
N = 100; 
TAU = 6; % {0, 1, 2, 3, 4, 5}
K = 95; % {0, 1, 2,..., 94}

% Generovani realizaci
Xk = zeros(N, M);
for i = 1:M
    for j = 2:N % X0 = 0
        Xk(j, i) = Xk(j - 1, i) + randn; 
    end
end

% Vykresleni 8 realizaci Wienerova procesu 
figure
plot(0:N-1, Xk(:, 1:8), 'LineWidth', 1);
xlabel('Cas [k]');
ylabel('X');
title('Vykresleni 8 realizaci Wienerova procesu ');

% Odhad parametru 
cov_Xk_odhad = zeros(TAU, K);
for k = 1:K
    for tau = 1:TAU
        Xk_x = Xk(k, :); % k-ty casovy okamzik
        Xk_y = Xk(k + tau - 1, :); % k + tau casovy okamzik
        EXk_x = mean(Xk_x); % stredni hodnota
        EXk_y = mean(Xk_y); % stredni hodnota
        for i = 1:length(Xk_x)
            cov_Xk_odhad(tau, k) = cov_Xk_odhad(tau, k) + ((Xk_x(i) - EXk_x) * (Xk_y(i) - EXk_y));
        end
        cov_Xk_odhad(tau, k) = cov_Xk_odhad(tau, k) / length(Xk_x);
    end
end

% Vypocet parametru 
cov_Xk_vypocet = zeros(TAU, K);
for tau = 1:TAU
    for k = 1:K
        cov_Xk_vypocet(tau, k) = k - 1; 
    end
end

% Vykresleni
figure
for j = 1:TAU
    hold on  
    plot(0:length(cov_Xk_odhad)-1, cov_Xk_odhad(j,:), 'LineWidth', 1.5); 
end

hold on
for u = 1:TAU
    hold on
    plot(0:length(cov_Xk_vypocet)-1, cov_Xk_vypocet(u, :), 'black:','LineWidth', 1.5);
end

Legend = cell(7,1);
Legend{1} = 'COV[X_k,X_k]'; 
Legend{2} = 'COV[X_k,X_{k+1}]';  
Legend{3} = 'COV[X_k,X_{k+2}]'; 
Legend{4} = 'COV[X_k,X_{k+3}]';
Legend{5} = 'COV[X_k,X_{k+4}]'; 
Legend{6} = 'COV[X_k,X_{k+5}]'; 
Legend{7} = 'Teoreticke hodnoty \newline jednotlivych COV[X_k,X_{k+\tau}]'; 

title('Autokovariancni funkce');
xlabel('Cas [k]');
legend(Legend, 'Location', 'southeast');

%% Priklad 3 
clc
clear all 
close all

% Zadane parametry 
M = 10^4;
N = 100; 
K = 100; % {0, 1, 2,..., 99}

% Generovani realizaci 
Xk = zeros(N, M);
Zk = zeros(N, M);
for i = 1:M
    Xk(1, i) = 1 + randn * sqrt(5);
    Zk(1, i) = 5 * Xk(1, i) + randn * sqrt(2);
    for j = 2:N
        Xk(j, i) = 0.95 * Xk(j - 1, i) + 0.5 * randn * sqrt(3);
        Zk(j, i) = 5 * Xk(j, i) + randn * sqrt(2);
    end
end

% Odhad parametru
E_Zk_odhad = zeros(1, K);
E_Xk_odhad = zeros(1, K);

var_Zk_odhad = zeros(1, K);
var_Xk_odhad = zeros(1, K);


for k = 1:K
    a = Xk(k, :);
    b = Zk(k, :);
    VARx = 0;
    VARz = 0;
    Ex = mean(a); 
    Ez = mean(b); 
    
    for i = 1:length(a)
        VARx = VARx + (a(i) - Ex)^2;
        VARz = VARz + (b(i) - Ez)^2;
    end
    
    VARx = VARx / length(a);
    VARz = VARz / length(b);
    
    E_Xk_odhad(k) = Ex;
    E_Zk_odhad(k) = Ez;
    
    var_Xk_odhad(k) = VARx;
    var_Zk_odhad(k) = VARz;
end

% Vypocet parametru
E_Xk_vypocet = zeros(1, K);
E_Zk_vypocet = zeros(1, K);

var_Xk_vypocet = zeros(1, K);
var_Zk_vypocet = zeros(1, K);

for k = 0:K
    E_Xk_vypocet(k+1) = 0.95^k;
    E_Zk_vypocet(k+1) = 5 * 0.95^k;
    
    suma = 0;
    for n = 0:k-1
        suma = suma + 0.95^(2*n);
    end
    
    var_Xk_vypocet(k+1) = 0.95^(2*k) * 5 + 0.5^2 * 3 * suma;
    var_Zk_vypocet(k+1) = var_Xk_vypocet(k+1) * 5^2 + 2;
end

% Vypocet ustalenych hodnot pro k->inf
k_ust = 100000; 

E_Xk_ust = 0.95^k_ust
E_Zk_ust = 5 * 0.95^k_ust
    
suma = 0;
for n = 0:k_ust
    suma = suma + 0.95^(2*n);
end
    
var_Xk_ust = 0.95^(2*k_ust) * 5 + 0.5^2 * 3 * suma
var_Zk_ust = var_Xk_ust * 5^2 + 2



% Vykresleni
% Stredni hodnoty X_k, Z_k
figure
hold on 
plot(0:length(E_Xk_odhad)-1, E_Xk_odhad, 'blue','LineWidth', 1.0);
plot(0:length(E_Xk_vypocet)-1, E_Xk_vypocet, 'black:','LineWidth', 1.5);
title('Stredni hodnoty X_k')
xlabel('Cas [k]')
ylabel('E[X_k]')
Legend = cell(2,1);
Legend{1} = 'Odhadovane stredni hodnoty X_k'; 
Legend{2} = 'Teoreticky vypoctene stredni hodnoty X_k';
legend(Legend); 

figure
hold on 
plot(0:length(E_Zk_odhad)-1, E_Zk_odhad, 'blue','LineWidth', 1.0);
plot(0:length(E_Zk_vypocet)-1, E_Zk_vypocet, 'black:','LineWidth', 1.5);
title('Stredni hodnoty Z_k')
xlabel('Cas [k]')
ylabel('E[Z_k]')
Legend = cell(2,1);
Legend{1} = 'Odhadovane stredni hodnoty Z_k'; 
Legend{2} = 'Teoreticky vypoctene stredni hodnoty Z_k';
legend(Legend); 

% Variance X_k, Z_k
figure
hold on 
plot(0:length(var_Xk_odhad)-1, var_Xk_odhad, 'blue','LineWidth', 1.0);
plot(0:length(var_Xk_vypocet)-1, var_Xk_vypocet, 'black:','LineWidth', 1.5);
title('Variance X_k')
xlabel('Cas [k]')
ylabel('VAR[X_k]')
Legend = cell(2,1);
Legend{1} = 'Odhadovane hodnoty varianci X_k'; 
Legend{2} = 'Teoreticky vypoctene hodnoty varianci X_k';
legend(Legend, 'Location','southeast'); 

figure
hold on 
plot(0:length(var_Zk_odhad)-1, var_Zk_odhad, 'blue','LineWidth', 1.0);
plot(0:length(var_Zk_vypocet)-1, var_Zk_vypocet, 'black:','LineWidth', 1.5);
title('Variance Z_k')
xlabel('Cas [k]')
ylabel('VAR[Z_k]')
Legend = cell(2,1);
Legend{1} = 'Odhadovane hodnoty varianci Z_k'; 
Legend{2} = 'Teoreticky vypoctene hodnoty varianci Z_k';
legend(Legend, 'Location','southeast'); 










