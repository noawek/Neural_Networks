%% introduction to neuronal networks - assignment 1:



% settimg some parameters:

Eps = 24;
R1 = 2;
R2 = linspace(0,22);



% resistors in series:

R_t_l = R1 + R2;
I_t_l = Eps./R_t_l;

figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(R2, I_t_l);
xlabel({'R_{2} [\Omega]' , '(resistence of the changing resistor in Ohm)'});
ylabel({'I_{T} [A]' , '(total current in the circuit in Ampere)'});
xlim([0 (max(R2) + 1)]);
xticks(0:1:(max(R2) + 1));
ylim([0 (max(I_t_l) + 1)]);
yticks(0:1:(max(I_t_l) + 1));
title({'', '\fontsize{13} \color{blue} - resistors are connected in series -', '', '\fontsize{14} \color{black} The total electric current flowing in an electronic circuit (I_{T})' , 'as a function of the changing resistence of one resistor (R_{2}) in that circuit',  ''});
hold off; 



% resistors in parallel:

R_t_p = 1./(1/R1 + 1./R2);
I_t_p = Eps./R_t_p;

figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(R2, I_t_p);
xlabel({'R_{2} [\Omega]' , '(resistence of the changing resistor in Ohm)'});
ylabel({'I_{T} [A]' , '(total current in the circuit in Ampere)'});
xlim([0 (max(R2) + 1)]);
xticks(0:1:(max(R2) + 1));
yticks(0:6:120);
asy = xline(0,'--r', 'LineWidth', 2, 'DisplayName','asymptote $(y \rightarrow \infty)$');
legend(asy, 'FontSize', 15, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title({'', '\fontsize{13} \color{blue} - resistors are connected in parallel -', '', '\fontsize{14} \color{black} The total electric current flowing in an electronic circuit (I_{T})' , 'as a function of the changing resistence of one resistor (R_{2}) in that circuit',  ''});
hold off;


