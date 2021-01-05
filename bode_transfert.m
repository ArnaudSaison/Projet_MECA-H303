function [] = bode_transfert(T, a, b, n, m_a, c_a, k_a, xi_a, beta)
    font_size = 11;
%     omega = logspace(a, b, n);
% 
%     G = 20*log10(abs(T(omega, m_a, c_a, k_a)));
%     legend_str = num2str(m_a, "m_a = %.2f") + ", " ...
%                + num2str(c_a, "c_a = %.2f") + ", " ...
%                + num2str(k_a, "k_a = %.2f");
%     semilogx(omega, G, 'DisplayName', char(legend_str)); 
%     grid on; hold on; legend('show');
%     xlabel("\omega (rad)"); ylabel("Mag (dB)");
    
    omega = linspace(a, b, n);

    G = 20*log10(T(omega, m_a, c_a, k_a));
    legend_str = num2str(m_a, "m_a = %.2f") + ", " ...
               + num2str(c_a, "c_a = %.2f") + ", " ...
               + num2str(k_a, "k_a = %.2f") + ", " ...
               + '\xi_a' + num2str(xi_a(m_a, c_a, k_a), " = %.2f") + ", " ...
               + '\beta' + num2str(beta(m_a, k_a), " = %.2f");
    plot(omega, G, 'DisplayName', legend_str);
    grid on; hold on; legend('show', 'Location','southwest', 'FontSize', font_size);
    xlabel("\omega (rad)", 'FontSize', font_size); ylabel("Mag (dB)", 'FontSize', font_size);  
end