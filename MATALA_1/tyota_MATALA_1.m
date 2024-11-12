%% introduction to neuronal networks - assignment 1:

clear;
clc;

% settimg some parameters:


% condutence of ions in [mS/cm^2]:
G_na = 120;
G_k = 36;
G_m = 0.3;  

% nernst potential of ions in [mV]:
E_na = 115;
E_k = -12;
E_m = 10.6;

% the membrane capacity in [(10^-6)*F*(cm^-2)]:
C_m = 1;

% time and deltaT in [mS]:
deltaT = 0.01;
Time = 50;
T = 0:deltaT:Time;

% Vrest and voltage in [mV]:
V_rest = 0;
V = zeros(size(T), 'like', T);
V(1) = V_rest;

% need fixing!!!!!




% the probability parameters n, m and h, and their rate constants (alpha and beta in [1/Sec]):

% n:
alpha_n = zeros(size(T), 'like', T);
beta_n = zeros(size(T), 'like', T);
n = zeros(size(T), 'like', T);
alpha_n(1) = (10-V(1))/(100*(exp((10-V(1))/10)-1));
beta_n(1) = 0.125*(exp((-V(1))/80));
n0 = alpha_n(1)/(alpha_n(1)+beta_n(1));
n(1) = n0;

% m:
alpha_m = zeros(size(T), 'like', T);
beta_m = zeros(size(T), 'like', T);
m = zeros(size(T), 'like', T);
alpha_m(1) = (25-V(1))/(10*(exp((25-V(1))/10)-1));
beta_m(1) = 4*(exp((-V(1))/18));
m0 = alpha_m(1)/(alpha_m(1)+beta_m(1));
m(1) = m0;

% h:
alpha_h = zeros(size(T), 'like', T);
beta_h = zeros(size(T), 'like', T);
h = zeros(size(T), 'like', T);
alpha_h(1) = 0.07*(exp((-V(1))/20));
beta_h(1) = 1/(exp((30-V(1))/10)+1);
h0 = alpha_h(1)/(alpha_h(1)+beta_h(1));
h(1) = h0;




% Iinj:

I_inj = zeros(size(T), 'like', T);
CurrentInj = 13.5; %*10^-6;   % [nA]
CurrentDuration_ms = 0.5;
CurrentStartTime_ms = 10;
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = CurrentStartTime+CurrentDuration;
I_inj(CurrentStartTime:CurrentEndTime) = CurrentInj;




dvdt = zeros(size(T), 'like', T);
for i = 1:(length(T)-1)
    
    % the probability parameters n, m and h, and their rate constants (alpha and beta):
    
    % n:
    alpha_n(i) = (10-V(i))/(100*(exp((10-V(i))/10)-1));
    beta_n(i) = 0.125*(exp((-V(i))/80));
    n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*deltaT;
    % m:
    alpha_m(i) = (25-V(i))/(10*(exp((25-V(i))/10)-1));
    beta_m(i) = 4*(exp((-V(i))/18));
    m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*deltaT;
    % h:
    alpha_h(i) = 0.07*(exp((-V(i))/20));
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*deltaT;
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;

end

V = V-65;
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(T, V, T, I_inj-90);
ylim([-100 70]);  % change to not manual!!
xlabel({'Time', '[mS]'});
ylabel({'Voltage', '[mV]'});
legend('FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
hold off;


figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(T, n, T, m, T, h);
xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
ylabel({'V', '[V]' , '(the voltage on the capacitor of the circuit in Volts)'});
legend('n', 'm', 'h', 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
hold off;


% %% introduction to neuronal networks - assignment 1:
% 
% 
% 
% %% settimg some parameters:
% 
% % condutence of ions in [mS/cm^2]:
% G_na = 120;
% G_k = 36;
% G_m = 0.3;  
% 
% % nernst potential of ions in [mV]:
% E_na = 115;
% E_k = -12;
% E_m = 115;
% 
% % the membrane capacity in [(10^-6)*F*(cm^-2)]:
% C_m = 1;
% 
% % time and deltaT in [mS]:
% deltaT = 0.01;
% Time = 100;
% T = 0:deltaT:Time;
% 
% % Vrest and voltage in [mV]:
% Vrest = 0;
% V = zeros(size(T), 'like', T);
% V(1) = Vrest;
% 
% % need fixing!!!!!
% 
% 
% %% the probability parameters n, m and h, and their rate constants (alpha and beta in [1/Sec]):
% 
% % n:
% alpha_n = zeros(size(T), 'like', T);
% beta_n = zeros(size(T), 'like', T);
% n = zeros(size(T), 'like', T);
% alpha_n(1) = (10-V(1))/(100*(exp((10-V(1))/10)-1));
% beta_n(1) = 0.125*(exp((-V(1))/80));
% n0 = alpha_n(1)/(alpha_n(1)+beta_n(1));
% n(1) = n0;
% 
% % m:
% alpha_m = zeros(size(T), 'like', T);
% beta_m = zeros(size(T), 'like', T);
% m = zeros(size(T), 'like', T);
% alpha_m(1) = (25-V(1))/(10*(exp((25-V(1))/10)-1));
% beta_m(1) = 4*(exp((-V(1))/18));
% m0 = alpha_m(1)/(alpha_m(1)+beta_m(1));
% m(1) = m0;
% 
% % h:
% alpha_h = zeros(size(T), 'like', T);
% beta_h = zeros(size(T), 'like', T);
% h = zeros(size(T), 'like', T);
% alpha_h(1) = 0.07*(exp((-V(1))/20));
% beta_h(1) = 1/(exp((30-V(1))/10)+1);
% h0 = alpha_h(1)/(alpha_h(1)+beta_h(1));
% h(1) = h0;
% 
% 
% %% the equation:
% 
% dvdt = (G_na*m^3*h*(E_na-V) + G_k*n^4*(E_k-V) + G_m*(V_rest-V) + I_inj(i))/C_m;
% 
% 
% %%
% 
% % Iinj:
% 
% Iinj = zeros(size(T), 'like', T);
% CurrentInj = 0.4; %*10^-6;   % [nA]
% CurrentDuration_ms = 0.5;
% CurrentStartTime_ms = 1;
% CurrentDuration = CurrentDuration_ms/deltaT;
% CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
% CurrentEndTime = CurrentStartTime+CurrentDuration;
% Iinj(CurrentStartTime:CurrentEndTime) = CurrentInj;
% 
% 
% %%
% 
% for i = 1:length(T)
%     
%     % the probability parameters n, m and h, and their rate constants (alpha and beta):
%     
%     % n:
%     alpha_n(i) = (10-V(i))/(100*(exp((10-V(i))/10)-1));
%     beta_n(i) = 0.125*(exp((-V(i))/80));
%     n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*deltaT;
%     % m:
%     alpha_m(i) = (25-V(i))/(10*(exp((25-V(i))/10)-1));
%     beta_m(i) = 4*(exp((-V(i))/18));
%     m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*deltaT;
%     % h:
%     alpha_h(i) = 0.07*(exp((-V(i))/20));
%     beta_h(i) = 1/(exp((30-V(i))/10)+1);
%     h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*deltaT;
%     dvdt(i) = (G_na*(m(i))^3*h(i)*(E_na-V(i)) + G_k*(n(i))^4*(E_k-V(i)) + G_m*(V_rest-V(i)) + I_inj(i))/C_m;
%     V(i+1) = V(i) + dvdt(i)*deltaT;
% 
% 
% 
% 
% end
% 
% 
% 
% %% the equation:
% 
% C_m*dvdt = G_na*m^3*h*(E_na-V) + G_k*n^4*(E_k-V) + G_m*(V_rest-V) + I_inj(i);
% 
% 
% V(i) = V(i-1) + dvdt(i)*deltaT;
% 
% 
% 
% R1 = 3000;
% R2 = 6000;
% C = 0.002;
% T_charge = linspace(0,20);
% T_discharge = linspace(0,30);
% 
% 
% % charging the capacitor:
% 
% I_c = (Eps/R1)*exp(-T_charge/(R1*C));
% V_c = Eps*C*(1-(exp(-T_charge/(R1*C))));
% 
% 
% % discharging the capacitor:
% 
% V1 = V_c(100);
% I_d = ((-V1)/R2)*exp(-T_discharge/(R2*C));
% V_d = V1*exp(-T_discharge/(R2*C));
% 
% 
% % current on the capacitor:
% 
% figure('Units','normalized','position',[0 0 1 1]);
% hold on;
% plot([T_charge, (T_discharge + 20)], [I_c, I_d]);
% xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
% ylabel({'I', '[A]' , '(the current on the capacitor of the circuit in Ampere)'});
% asy_current = xline(20,'--r', 'LineWidth', 0.5);
% yline(0,'--k', 'LineWidth', 0.5);
% legend(asy_current, '$charging \rightarrow discharging$', 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
% title({'\fontsize{14} \color{blue} - the current on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (current VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
% hold off; 
% 
% 
% % voltage on the capacitor:
% 
% figure('Units','normalized','position',[0 0 1 1]);
% hold on;
% plot([T_charge, (T_discharge + 20)], [V_c, V_d]);
% xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
% ylabel({'V', '[V]' , '(the voltage on the capacitor of the circuit in Volts)'});
% asy_voltage = xline(20,'--r', 'LineWidth', 0.5, 'DisplayName','$charging \rightarrow discharging$');
% legend(asy_voltage, 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
% title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
% hold off;

