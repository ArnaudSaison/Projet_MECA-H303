clear all; 
%close all;
figure;
clc;

for i = [0.1, 1, 10]
% Données
k = 30;
m = 1;
c = i;

omega_n = sqrt(k/m)
xi = c / (2 * m * omega_n)
omega_d = omega_n * sqrt(1 - xi^2)

% Diagramme de Bode : réponse fréquentielle
H =@(omega) (1./(-m.*omega.^2 + 1i.*c.*omega + k));

% Réponse impulsionnelle
x =@(t)  exp(-xi.*omega_n.*t).*sin(omega_d.*t)/(m.*omega_d);

% Traçage
dia_bode(H, x, 1000, 0, 1.5, 50, num2str(c, "c = %.2f"));
end