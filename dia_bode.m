function [] = dia_bode(h, y, n, a, bh, by, legend_str)
    % ---------------------------------------------------------------------    
    % Tracé des courbes de Bode et réponse impulsionnelle
    % ---------------------------------------------------------------------    
    % h = transmittance isochrone
    % y = réponse impulsionnelle
    % n = nombre de pas pour les graphiques
    % a = 0 (borne inférieure sur les graphiques
    % bh = borne supérieure en 10^b pour les diagrammes de Bode
    % by = borne sup en t pour le graphe de la réponse impulsionnelle
    % legend_str = légende qui apparaît pour tous les graphes
    %
    % Pour afficher plusieurs tracés sur ces grpahes, supprimer toute ligne
    % "close all" et simplement appeler plusieurs fois la fonction (peut
    % être fait en plusieurs scripts)
    % ---------------------------------------------------------------------
    
    % Axe des fréquences (pulsations en microseconde^{-1})
    omega = logspace(a, bh, n);
    H = h(omega);
    
    % module de H
    G = 20*log10(abs(H));

    % Argument (principal) de H
    Phi = angle(H);

    % Mise en graphique
    % Bode : amplitude
    subplot(3,1,1); semilogx(omega, G, 'DisplayName', char(legend_str)); 
    grid on; hold on; legend('show');
    xlabel("\omega (rad)"); ylabel("Mag (dB)");

    % Bode : phase
    subplot(3,1,2); semilogx(omega, Phi * 180/pi, 'DisplayName', legend_str); 
    grid on; hold on; legend('show');
    xlabel("\omega (rad)"); ylabel("Phase (deg)");
    
    % Réponse impulsionnelle
    t = linspace(a, by, n);
    Y = y(t);
    subplot(3,1,3); plot(t, Y, 'DisplayName', legend_str); 
    grid on; hold on; legend('show');
    xlabel("t (s)"); ylabel("Amplitude (rad)");
end