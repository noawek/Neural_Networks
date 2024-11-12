%% Neural Networks - assignment 2:



% settimg some parameters:

Eps = 30;
R1 = 3000;
R2 = 6000;
C = 0.002;
T_charge = linspace(0,20);
T_discharge = linspace(0,30);


% charging the capacitor:

I_c = (Eps/R1)*exp(-T_charge/(R1*C));
V_c = Eps*C*(1-(exp(-T_charge/(R1*C))));


% discharging the capacitor:

V1 = V_c(100);
I_d = ((-V1)/R2)*exp(-T_discharge/(R2*C));
V_d = V1*exp(-T_discharge/(R2*C));


% current on the capacitor:

figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot([T_charge, (T_discharge + 20)], [I_c, I_d]);
xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
ylabel({'I', '[A]' , '(the current on the capacitor of the circuit in Ampere)'});
asy_current = xline(20,'--r', 'LineWidth', 0.5);
yline(0,'--k', 'LineWidth', 0.5);
legend(asy_current, '$charging \rightarrow discharging$', 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title({'\fontsize{14} \color{blue} - the current on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (current VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
hold off; 


% voltage on the capacitor:

figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot([T_charge, (T_discharge + 20)], [V_c, V_d]);
xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
ylabel({'V', '[V]' , '(the voltage on the capacitor of the circuit in Volts)'});
asy_voltage = xline(20,'--r', 'LineWidth', 0.5, 'DisplayName','$charging \rightarrow discharging$');
legend(asy_voltage, 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
hold off;



