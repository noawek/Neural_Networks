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

alpha_n(1) = (10-V(1))/(100*(exp((10-V(1))/10)-1));
beta_n(1) = 0.125*(exp((-V(1))/80));
n0 = alpha_n(1)/(alpha_n(1)+beta_n(1));
n = n0;

% m:

alpha_m(1) = (25-V(1))/(10*(exp((25-V(1))/10)-1));
beta_m(1) = 4*(exp((-V(1))/18));
m0 = alpha_m(1)/(alpha_m(1)+beta_m(1));
m = m0;

% h:

alpha_h(1) = 0.07*(exp((-V(1))/20));
beta_h(1) = 1/(exp((30-V(1))/10)+1);
h0 = alpha_h(1)/(alpha_h(1)+beta_h(1));
h = h0;




% Iinj:


CurrentInj = flip(0:70); %*10^-6;   % [nA]
CurrentStartTime_ms = 10;
CurrentDuration_ms = 0.1:0.1:(Time-CurrentStartTime_ms);    %% have to fix!!!
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = round(CurrentStartTime+CurrentDuration);
I_inj = zeros(size(T), 'like', T);
I_th = zeros(1, length(CurrentDuration));

%%
for k = 1:length(CurrentDuration)
    for j = 1:length(CurrentInj)
        I_inj(CurrentStartTime:CurrentEndTime(k)) = CurrentInj(j);
        dvdt = zeros(size(T), 'like', T);
        for i = 1:(length(I_inj)-1)
            % the probability parameters n, m and h, and their rate constants (alpha and beta):
            dvdt = (G_na*((m)^3)*h*(E_na-V(i)) + G_k*((n)^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;

            % n:
            alpha_n = (10-V(i))/(100*(exp((10-V(i))/10)-1));
            beta_n = 0.125*(exp((-V(i))/80));
            n = n + (alpha_n*(1-n)-beta_n*n)*deltaT;
            % m:
            alpha_m = (25-V(i))/(10*(exp((25-V(i))/10)-1));
            beta_m = 4*(exp((-V(i))/18));
            m = m + (alpha_m*(1-m)-beta_m*m)*deltaT;
            % h:
            alpha_h = 0.07*(exp((-V(i))/20));
            beta_h = 1/(exp((30-V(i))/10)+1);
            h = h + (alpha_h*(1-h)-beta_h*h)*deltaT;
            
            V(i+1) = V(i) + dvdt*deltaT;
        end
        TH = 30;
        Above = V>=TH;
        PlacesAbove = find(Above);
        TheCross = V(PlacesAbove-1)<TH;
        TimesCrossing = length(PlacesAbove(TheCross)); 
        if TimesCrossing == 0 && I_th(k) == 0
            I_th(k) = CurrentInj(j-1);
        end
        n(1) = n0;
        m(1) = m0;
        h(1) = h0;
        V(1) = V_rest;
    end
end

% %% introduction to neuronal networks - assignment 1:
% 
% clear;
% clc;
% 
% % settimg some parameters:
% 
% 
% % condutence of ions in [mS/cm^2]:
% G_na = 120;
% G_k = 36;
% G_m = 0.3;  
% 
% % nernst potential of ions in [mV]:
% E_na = 115;
% E_k = -12;
% E_m = 10.6;
% 
% % the membrane capacity in [(10^-6)*F*(cm^-2)]:
% C_m = 1;
% 
% % time and deltaT in [mS]:
% deltaT = 0.01;
% Time = 50;
% T = 0:deltaT:Time;
% 
% % Vrest and voltage in [mV]:
% V_rest = 0;
% V = zeros(size(T), 'like', T);
% V(1) = V_rest;
% 
% % need fixing!!!!!
% 
% 
% 
% 
% % the probability parameters n, m and h, and their rate constants (alpha and beta in [1/Sec]):
% 
% % n:
% n = zeros(size(T), 'like', T);
% alpha_n(1) = (10-V(1))/(100*(exp((10-V(1))/10)-1));
% beta_n(1) = 0.125*(exp((-V(1))/80));
% n0 = alpha_n(1)/(alpha_n(1)+beta_n(1));
% n(1) = n0;
% 
% % m:
% m = zeros(size(T), 'like', T);
% alpha_m(1) = (25-V(1))/(10*(exp((25-V(1))/10)-1));
% beta_m(1) = 4*(exp((-V(1))/18));
% m0 = alpha_m(1)/(alpha_m(1)+beta_m(1));
% m(1) = m0;
% 
% % h:
% h = zeros(size(T), 'like', T);
% alpha_h(1) = 0.07*(exp((-V(1))/20));
% beta_h(1) = 1/(exp((30-V(1))/10)+1);
% h0 = alpha_h(1)/(alpha_h(1)+beta_h(1));
% h(1) = h0;
% 
% 
% 
% 
% % Iinj:
% 
% 
% CurrentInj = flip(0:70); %*10^-6;   % [nA]
% CurrentStartTime_ms = 10;
% CurrentDuration_ms = 0.1:0.1:(Time-CurrentStartTime_ms);    %% have to fix!!!
% CurrentDuration = CurrentDuration_ms/deltaT;
% CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
% CurrentEndTime = round(CurrentStartTime+CurrentDuration);
% I_inj = zeros(size(T), 'like', T);
% I_th = zeros(1, length(CurrentDuration));
% 
% %%
% for k = 1:length(CurrentDuration)
%     for j = 1:length(CurrentInj)
%         I_inj(CurrentStartTime:CurrentEndTime(k)) = CurrentInj(j);
%         dvdt = zeros(size(T), 'like', T);
%         for i = 1:(length(I_inj)-1)
%             % the probability parameters n, m and h, and their rate constants (alpha and beta):
%             
%             % n:
%             n(i+1) = n(i) + (((10-V(i))/(100*(exp((10-V(i))/10)-1)))*(1-n(i))-(0.125*(exp((-V(i))/80)))*n(i))*deltaT;
%             % m:
%             m(i+1) = m(i) + (((25-V(i))/(10*(exp((25-V(i))/10)-1)))*(1-m(i))-(4*(exp((-V(i))/18)))*m(i))*deltaT;
%             % h:
%             h(i+1) = h(i) + ((0.07*(exp((-V(i))/20)))*(1-h(i))-(1/(exp((30-V(i))/10)+1))*h(i))*deltaT;
%             
%             V(i+1) = V(i) + ((G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m)*deltaT;
%         end
%         TH = 30;
%         Above = V>=TH;
%         PlacesAbove = find(Above);
%         TheCross = V(PlacesAbove-1)<TH;
%         TimesCrossing = length(PlacesAbove(TheCross)); 
%         if TimesCrossing == 0 && I_th(k) == 0
%             I_th(k) = CurrentInj(j-1);
%         end
%         n(1) = n0;
%         m(1) = m0;
%         h(1) = h0;
%         V(1) = V_rest;
%     end
% end
% 
% 
% 
% 
% 
% % %% introduction to neuronal networks - assignment 1:
% % 
% % clear;
% % clc;
% % 
% % % settimg some parameters:
% % 
% % 
% % % condutence of ions in [mS/cm^2]:
% % G_na = 120;
% % G_k = 36;
% % G_m = 0.3;  
% % 
% % % nernst potential of ions in [mV]:
% % E_na = 115;
% % E_k = -12;
% % E_m = 10.6;
% % 
% % % the membrane capacity in [(10^-6)*F*(cm^-2)]:
% % C_m = 1;
% % 
% % % time and deltaT in [mS]:
% % deltaT = 0.01;
% % Time = 50;
% % T = 0:deltaT:Time;
% % 
% % % Vrest and voltage in [mV]:
% % V_rest = 0;
% % V = zeros(size(T), 'like', T);
% % V(1) = V_rest;
% % 
% % % need fixing!!!!!
% % 
% % 
% % 
% % 
% % % the probability parameters n, m and h, and their rate constants (alpha and beta in [1/Sec]):
% % 
% % % n:
% % alpha_n = zeros(size(T), 'like', T);
% % beta_n = zeros(size(T), 'like', T);
% % n = zeros(size(T), 'like', T);
% % alpha_n(1) = (10-V(1))/(100*(exp((10-V(1))/10)-1));
% % beta_n(1) = 0.125*(exp((-V(1))/80));
% % n0 = alpha_n(1)/(alpha_n(1)+beta_n(1));
% % n(1) = n0;
% % 
% % % m:
% % alpha_m = zeros(size(T), 'like', T);
% % beta_m = zeros(size(T), 'like', T);
% % m = zeros(size(T), 'like', T);
% % alpha_m(1) = (25-V(1))/(10*(exp((25-V(1))/10)-1));
% % beta_m(1) = 4*(exp((-V(1))/18));
% % m0 = alpha_m(1)/(alpha_m(1)+beta_m(1));
% % m(1) = m0;
% % 
% % % h:
% % alpha_h = zeros(size(T), 'like', T);
% % beta_h = zeros(size(T), 'like', T);
% % h = zeros(size(T), 'like', T);
% % alpha_h(1) = 0.07*(exp((-V(1))/20));
% % beta_h(1) = 1/(exp((30-V(1))/10)+1);
% % h0 = alpha_h(1)/(alpha_h(1)+beta_h(1));
% % h(1) = h0;
% % 
% % 
% % 
% % 
% % % Iinj:
% % 
% % CurrentInj = flip(0:30); %*10^-6;   % [nA]
% % CurrentStartTime_ms = 10;
% % CurrentDuration_ms = 0:deltaT:(Time-CurrentStartTime_ms);
% % CurrentDuration = CurrentDuration_ms/deltaT;
% % CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
% % CurrentEndTime = CurrentStartTime+CurrentDuration;
% % I_inj = zeros(length(CurrentDuration), length(T), length(CurrentInj));
% % 
% % for k = 1:length(CurrentDuration)
% %     for j = 1:length(CurrentInj)
% %         I_inj(k, CurrentStartTime:CurrentEndTime(k), j) = CurrentInj(j);
% %         dvdt = zeros(size(T), 'like', T);
% %         for i = 1:(length(T)-1)
% %             % the probability parameters n, m and h, and their rate constants (alpha and beta):
% %             
% %             % n:
% %             alpha_n(i) = (10-V(i))/(100*(exp((10-V(i))/10)-1));
% %             beta_n(i) = 0.125*(exp((-V(i))/80));
% %             n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*deltaT;
% %             % m:
% %             alpha_m(i) = (25-V(i))/(10*(exp((25-V(i))/10)-1));
% %             beta_m(i) = 4*(exp((-V(i))/18));
% %             m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*deltaT;
% %             % h:
% %             alpha_h(i) = 0.07*(exp((-V(i))/20));
% %             beta_h(i) = 1/(exp((30-V(i))/10)+1);
% %             h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*deltaT;
% %             
% %             dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
% %             V(i+1) = V(i) + dvdt(i)*deltaT;
% %         end
% %         TH = 6.85;
% %         Above = V>=TH;
% %         PlacesAbove = find(Above);
% %         TheCross = V(PlacesAbove-1)<TH;
% %         TimesCrossing = length(PlacesAbove(TheCross)); 
% %         if TimesCrossing == 0
% %             I_th(k) = CurrentInj(j-1);
% %         end
% %         
% %     end
% % end
% % 
% % 
% % 
% % dvdt = zeros(size(T), 'like', T);
% % for i = 1:(length(T)-1)
% %     
% %     % the probability parameters n, m and h, and their rate constants (alpha and beta):
% %     
% %     % n:
% %     alpha_n(i) = (10-V(i))/(100*(exp((10-V(i))/10)-1));
% %     beta_n(i) = 0.125*(exp((-V(i))/80));
% %     n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*deltaT;
% %     % m:
% %     alpha_m(i) = (25-V(i))/(10*(exp((25-V(i))/10)-1));
% %     beta_m(i) = 4*(exp((-V(i))/18));
% %     m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*deltaT;
% %     % h:
% %     alpha_h(i) = 0.07*(exp((-V(i))/20));
% %     beta_h(i) = 1/(exp((30-V(i))/10)+1);
% %     h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*deltaT;
% %     dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
% %     V(i+1) = V(i) + dvdt(i)*deltaT;
% % 
% % end
% % 
% % V = V-65;
% % figure('Units','normalized','position',[0 0 1 1]);
% % hold on;
% % plot(T, V, T, I_inj-90);
% % ylim([-100 70]);  % change to not manual!!
% % xlabel({'Time', '[mS]'});
% % ylabel({'Voltage', '[mV]'});
% % legend('FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
% % title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
% % hold off;
% % 
% % 
% % figure('Units','normalized','position',[0 0 1 1]);
% % hold on;
% % plot(T, n, T, m, T, h);
% % xlabel({'Time', '[seconds]' , '(20 seconds of charging the capacitor, following by 30 seconds of discharging it)'});
% % ylabel({'V', '[V]' , '(the voltage on the capacitor of the circuit in Volts)'});
% % legend('n', 'm', 'h', 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
% % title({'\fontsize{14} \color{blue} - the voltage on a capacitor while charging it for 20 seconds, then discharging for 30 seconds -', '\fontsize{12} (voltage VS time)', '\fontsize{10} \color{black} R_{charge} =  3000 \Omega  ,        R_{discharge} =  6000 \Omega  ,        C  =  2 mF  ,        \fontsize{18} \epsilon \fontsize{10}  =  30 V  ,       V1_{t = 20}  =  0.0579 V'});
% % hold off;
% 
