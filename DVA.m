clear all; close all; font_size = 11; subplots = true;

%% Constantes
g = 9.81;

%% Données du problème
m = 1;
c = 0.1;
k = 730;

% Paramètre à donner pour résoudre les équations de Dan Hartog
% Les paramètres que l'on cherche sont c_a et k_a
lambda_max = 0.1;
m_a = lambda_max * m;

% if subplots == true; subplot(2,3,1); else; figure; end
% lambda_values = linspace(0,lambda_max, 1000);
% H_max =@(lambda) (1 + lambda) * sqrt((2 + lambda) / lambda);
% plot(lambda_values, H_max(lambda_values), '-r');
% title("Maximum", 'FontSize', font_size);

%% Paramètres des représentations graphiques
a = 15;
b = 35;
n = 1000;

%% Paramètres réduits et fonctions
% indépendants des paramètres du DVA
omega_p = sqrt(k / m);               % fréquence propre du système initial
xi_p = c / (2 * m * omega_p);        % coefficient d'amortissement initial
r =@(omega) omega / omega_p;         % rapport fréq d'excitation / système

% dépendants des paramètres du DVA
lambda =@(m_a) m_a / m;                                      % rapport de masses
omega_a =@(m_a, k_a) abs(sqrt(k_a / m_a));                   % fréquence propre du DVA
xi_a =@(m_a, c_a, k_a) c_a / (2 * m_a * omega_a(m_a, k_a));  % coefficient d'amortissement du DVA

% combinaisons de paramètres du DVA
beta =@(m_a, k_a) omega_a(m_a, k_a) / omega_p;     % rapport fréq DVA / système
eta =@(m_a, c_a, k_a) xi_p / xi_a(m_a, c_a, k_a);  % rapport des coefficients d'amortissement

% Fonction de transfert en fonction de omega et des paramètres
T =@(omega, m_a, c_a, k_a) sqrt(( (2 .* xi_a(m_a, c_a, k_a) .* beta(m_a, k_a) .* r(omega)).^2 + (beta(m_a, k_a).^2 - r(omega).^2).^2 )./( (r(omega).^2 .* (eta(m_a, c_a, k_a) + beta(m_a, k_a) + lambda(m_a) .* beta(m_a, k_a)) - beta(m_a, k_a) .* (1 + eta(m_a, c_a, k_a) .* beta(m_a, k_a))).^2 .* (2 .* xi_a(m_a, c_a, k_a) .* r(omega)).^2 + ((1 - r(omega).^2) .* (beta(m_a, k_a).^2 - r(omega).^2) - beta(m_a, k_a) .* r(omega).^2 .* (lambda(m_a) .* beta(m_a, k_a) + 4 .* eta(m_a, c_a, k_a) .* xi_a(m_a, c_a, k_a).^2)).^2));

%% Equal peak procedure
% Les deux équations de Den Hartog forment un système avec deux inconnues
F1 =@(m_a, c_a, k_a) beta(m_a, k_a) - 1 / (1 + lambda(m_a));
F2 =@(m_a, c_a, k_a) xi_a(m_a, c_a, k_a)^2 - 3 * lambda(m_a) / (8 * (1 + lambda(m_a)));

F =@(m_a, c_a, k_a) ([F1(m_a, c_a, k_a); F2(m_a, c_a, k_a)]);

% Vecteur servant de point de départ à la résolution numérique
x0 = [2; 100];

% Résolution
sol = fsolve(@(x) F(m_a, x(1), x(2)), x0, optimoptions('fsolve','Display','iter'));

c_a = sol(1);
k_a = sol(2);

%% diagramme de bode de la fonction de transfert
if subplots == true; subplot(2,3,3); else; figure; end

bode_transfert(T, 0, 50, n, m_a, c_a, k_a, xi_a, beta);
title("Diagramme de Bode de la solution", 'FontSize', font_size);
set(findall(gca, 'Type', 'Line'),'LineWidth',2)

%% Comparaison aux cas quasi-nul et sans DVA
hold off;
if subplots == true; subplot(2,3,5); else; figure; end

% On prend une valeur très petite de c_a pour simuler un amortissement nul
m_a_nul = 0.0001;
c_a_nul = 0.0001;
k_a_nul = 0.0001;

c_a_faible = 0.001;

bode_transfert(T, 0, 50, n, m_a, c_a, k_a, xi_a, beta);
title("Comparaison avec les cas sans DVA et à amortissement très faible", 'FontSize', font_size);
set(findall(gca, 'Type', 'Line'),'LineWidth',2)

% pas de DVA
hold on;
bode_transfert(T, 0, 50, n, m_a_nul, c_a_nul, k_a_nul, xi_a, beta);

% très faible amortissement
hold on;
bode_transfert(T, 0, 50, n, m_a, c_a_faible, k_a, xi_a, beta);

% légende des plots
legend("DVA correctement dimensionné", "sans DVA", "très faible amortissement");


%% Représentation d'autres valeurs pour comparaison
% Range dans lequel on fait varier les paramètres
range = 0.15;

% Calcul de P et Q
[P, Q] = calc_P_Q(T, m_a, c_a, k_a, range, b, n);

% -------------------------------------------------------------------------
% ETAPE 1 : Valeurs dans range autour des bonnes valeurs
% -------------------------------------------------------------------------
% nouvelle figure
hold off;
if subplots == true; subplot(2,3,1); else; figure; end

% On retrace le graphe de la soltuion
bode_transfert(T, a, b, n, m_a, c_a, k_a, xi_a, beta);
title("Etape 1 : P et Q à la même hauteur", 'FontSize', font_size);
set(findall(gca, 'Type', 'Line'),'LineWidth',2)

% cas 1
bode_transfert(T, a, b, n, m_a * (1 + range), c_a, k_a, xi_a, beta);

% cas 2
bode_transfert(T, a, b, n, m_a * (1 - range), c_a, k_a, xi_a, beta);

% P et Q
xline(P(1), '--k', {'\omega_P'},'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', 'FontSize', font_size);
xline(Q(1), '--k', {'\omega_Q'},'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', 'FontSize', font_size);

% -------------------------------------------------------------------------
% ETAPE 2 : Valeurs dans range autour des bonnes valeurs
% -------------------------------------------------------------------------
% nouvelle figure
hold off;
if subplots == true; subplot(2,3,2); else; figure; end

% On retrace le graphe de la soltuion
bode_transfert(T, a, b, n, m_a, c_a, k_a, xi_a, beta);
title("Etape 2 : P et Q maximums", 'FontSize', font_size);
set(findall(gca, 'Type', 'Line'),'LineWidth',2)

% cas 1
bode_transfert(T, a, b, n, m_a, c_a * (1 + range), k_a, xi_a, beta);

% cas 2
bode_transfert(T, a, b, n, m_a, c_a * (1 - range), k_a, xi_a, beta);

% cas 3
bode_transfert(T, a, b, n, m_a, c_a * (1 + 2 * range), k_a, xi_a, beta);

% cas 4
bode_transfert(T, a, b, n, m_a, c_a * (1 - 2 * range), k_a, xi_a, beta);

% P et Q
xline(P(1), '--k', {'\omega_P'},'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', 'FontSize', font_size);
xline(Q(1), '--k', {'\omega_Q'},'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', 'FontSize', font_size);


%% Analogie des paramètres
c_p = c_a;
m_p = m_a;
l = m_p * g / k_a;

% Affichage des résultats
disp("Paramètres pour le système équivalent avec deux masses :" ...
     + newline + "m_a = " + num2str(m_a, "%.2f") ...
     + newline + "c_a = " + num2str(c_a, "%.2f") ...
     + newline + 'k_a = ' + num2str(k_a, "%.4f"));
disp(newline + "Paramètres pour le système avec pendule :" ...
     + newline + "m_p = " + num2str(m_p, "%.2f") ...
     + newline + "c_p = " + num2str(c_p, "%.2f") ...
     + newline + "l = " + num2str(l, "%.5f"));
 
%% Comparaison avec d'autres valeurs de rapport de masses
hold off;
if subplots == true; subplot(2,3,4); else; figure; end

% valeurs du rapport de masses
lambdas = [0.1 0.08 0.06 0.04 0.02];

for lambda_max = lambdas
    % Conversion en masse
    m_a = lambda_max * m;
    
    % résolution pour ce rapport de masses
    sol = fsolve(@(x) F(m_a, x(1), x(2)), x0, optimset("Display", "off"));
    c_a = sol(1);
    k_a = sol(2);
    
    omega = linspace(a, b, n);
    G = 20*log10(T(omega, m_a, c_a, k_a));
    legend_str = "\lambda = " + num2str(lambda_max, "%.2f");
    plot(omega, G, 'DisplayName', legend_str);
    hold on;
end

grid on; legend('show', 'Location','southwest', 'FontSize', font_size);
xlabel("\omega (rad)", 'FontSize', font_size); ylabel("Mag (dB)", 'FontSize', font_size); 
title("Comparaison des \lambda", 'FontSize', font_size);
