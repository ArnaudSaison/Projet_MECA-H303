clear all; 
close all; 
clc;

for i = [0.01, 0.5, 1]
% Données
L = 0.2;
m = 1;
c_p = i;
g = 9.81;

omega_n = sqrt(g/L);
xi = c_p / (2 * m * omega_n);
omega_d = omega_n * sqrt(1 - xi^2);

% Diagramme de Bode : réponse fréquentielle
H =@(omega) 1./(-m.*L^2.*omega + 1i.*c_p.*L.^2.*omega + m.*g.*L);

% Réponse impulsionnelle
theta =@(t) exp(-xi.*omega_n.*t) /(m*L^2*omega_d) .* sin(omega_d.*t);

% Traçage
dia_bode(H, theta, 1000, 0, 3, 30, num2str(c_p, "c_p = %.2f"));
end