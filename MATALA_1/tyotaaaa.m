


%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%



clear;
clc;




%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
Time = 100;
T = 0:deltaT:Time;


% V_rest and voltage in [mV]:
V_rest = 0;
V = zeros(size(T), 'like', T);
V(1) = V_rest;
dvdt = zeros(size(T), 'like', T);


% I_inj (an external current injected to the nueron) in [nA]:         
I_inj = zeros(size(T), 'like', T);
CurrentInj = 15; 
CurrentDuration_ms = 0.5;
CurrentStartTime_ms = 10;
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = CurrentStartTime+CurrentDuration;
I_inj(CurrentStartTime:CurrentEndTime) = CurrentInj;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n, m & h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the probability parameters: n, m & h, and their voltage dependent rate constants (alpha and beta in [1/Sec]):


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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% voltage step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:(length(T)-1)
    
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
    % voltage:
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;

end


% a figure showing the existence of an action potential after a short pulse: 
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(4,1,1:3);
plot(T, V);
title('Voltage VS Time');
ylim([min(V)-10 max(V)+10]);                           
xlabel('Time [mS]');
ylabel({'Voltage', '[mV]'});
subplot(4,1,4);
plot(T, I_inj, Color = 'red');
title({'','Injected Current VS Time'});
xlabel('Time [mS]');
ylabel({'Current', 'I', '[nA]'});
ylim([min(I_inj)-10 max(I_inj)+10]);                           
sgtitle({'            an Action Potential due to a Short Pulse' , ''})
hold off;


% figure of the dynamics of the model's parameters:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(4,1,1:3);
plot(T, n, T, m, T, h);
xlabel('Time [mS]');
ylabel('Gating Variables');
legend('n', 'm', 'h', 'FontSize', 12, 'Location','best', 'Interpreter','latex', 'Box', 'off');
title('\fontsize{11} The dynamics of \fontsize{17} m, n & h \fontsize{11} over time');
subplot(4,1,4);
plot(T, I_inj, Color = 'green');
title({'','Injected Current VS Time'});
xlabel('Time [mS]');
ylabel({'Current', 'I', '[nA]'});
ylim([min(I_inj)-10 max(I_inj)+10]); 
sgtitle({'            Dynamics of \fontsize{20} m, n & h \fontsize{14} due to a Short Pulse' , ''})
hold off;




%%   QUESTION 2   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters again %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% I_inj (a vector of external currents varying in their intensities) in [nA]:         
CurrentInj = flip(linspace(0, 110, 220)); 
CurrentStartTime_ms = 10;
CurrentDuration_ms = linspace(0.05, (Time-CurrentStartTime_ms));   
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = round(CurrentStartTime+CurrentDuration);
I_inj = zeros(size(T), 'like', T);
I_th = zeros(1, length(CurrentDuration));



for k = 1:length(CurrentDuration)
    for j = 1:length(CurrentInj)
        I_inj(CurrentStartTime:CurrentEndTime(k)) = CurrentInj(j);
        dvdt = zeros(size(T), 'like', T);
        n(1) = n0;
        m(1) = m0;
        h(1) = h0;
        V(1) = V_rest;
        for i = 1:(length(I_inj)-1)
            
            % the gating variables: n, m and h, and their rate constants (alpha and beta):
            
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
        TH = 30;
        Above = V>=TH;
        PlacesAbove = find(Above);
        TheCross = V(PlacesAbove-1)<TH;
        TimesCrossing = length(PlacesAbove(TheCross)); 
        if TimesCrossing == 0 && I_th(k) == 0
            I_th(k) = CurrentInj(j-1);
        end
    end
end



figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(CurrentDuration_ms, I_th);
title('\fontsize{14} The threshold current of the system for pulses in different lengths');
xlabel('Pulse Duration [mS]');
ylabel({'Threshold Current', 'I', '[nA]'});




%%   QUESTION 3   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters again %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% time and deltaT in [mS]:
deltaT = 0.01;
Time = 2000;
T = 0:deltaT:Time;


% I_inj (an external current injected to the nueron) in [nA]:     
CurrentInj = 0.5:0.5:100; 
I_inj = zeros(size(T), 'like', T);
TimesCrossing = zeros(size(CurrentInj), 'like', CurrentInj);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% voltage step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for j = 1:length(CurrentInj)
    
    I_inj(:) = CurrentInj(j);
    n(1) = n0;
    m(1) = m0;
    h(1) = h0;
    V(1) = V_rest;
    for i = 1:(length(T)-1)
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
    % voltage:
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;
    end
    TH = 30;
    Above = V>=TH;
    PlacesAbove = find(Above);
    TheCross = V(PlacesAbove-1)<TH;
    TimesCrossing(j) = length(PlacesAbove(TheCross)); 

end



figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(CurrentInj,TimesCrossing/2);
title('\fontsize{14} The relationship between the intensity of a costant injected current and the firing rate of the neuron');
xlabel({'Current Intensity', 'I', '[nA]'});
ylabel({'Firing Rate', '[Hz]'});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% I_inj (an external current injected to the nueron) in [nA]:          fix!!
CurrentInj = [4, 25, 75, 100];  %*10^-6;   % [nA]
I_inj = zeros(length(CurrentInj), length(T));
TimesCrossing = zeros(size(CurrentInj), 'like', CurrentInj);


figure('Units','normalized','position',[0 0 1 1]); 
sgtitle({'Examples for firing rates (shown as voltage VS time)', 'due to different intensities of the constant injected current:', ''})
for j = 1:length(CurrentInj)
    
    I_inj(:) = CurrentInj(j);
    n(1) = n0;
    m(1) = m0;
    h(1) = h0;
    V(1) = V_rest;
    
    for i = 1:(length(T)-1)
    
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
    % voltage:
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;
    
    end

    TH = 30;
    Above = V>=TH;
    PlacesAbove = find(Above);
    TheCross = V(PlacesAbove-1)<TH;
    TimesCrossing(j) = length(PlacesAbove(TheCross)); 

    % a figure showing the existence of an action potential after a short pulse: 
    hold on;
    subplot(2,2,j)
    plot(T(1:10001), V(1:10001));
    xlabel({'Time [mS]'});
    ylabel({'Voltage', '[mV]'});
    title(['Firing Rate of: ', num2str(round(TimesCrossing(j)/2)), ' Hz   when the intensity of the current is: ', num2str(CurrentInj(j)), ' nA']);
    hold on;

end
hold off;




%%   QUESTION 4   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% time and deltaT in [mS]:
deltaT = 0.01;
Time = 100;
T = 0:deltaT:Time;


% V_rest and voltage in [mV]:
V = zeros(size(T), 'like', T);
V(1) = V_rest;
dvdt = zeros(size(T), 'like', T);


% DEFINING REFRACTORY PERIOD (time after stimulus in mS)
REF_ABS_ms = 1;
REF_REL_ms = 20;
REF_ABS = REF_ABS_ms/deltaT;
REF_REL = REF_REL_ms/deltaT;


% I_inj (an external current injected to the nueron) in [nA]:         
I_inj = zeros(size(T), 'like', T);
CurrentInj = 15; 
CurrentDuration_ms = 0.5;
CurrentStartTime_ms = 10;
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = CurrentStartTime+CurrentDuration;
I_inj([CurrentStartTime:CurrentEndTime, CurrentStartTime+REF_REL:CurrentEndTime+REF_REL]) = CurrentInj;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n, m & h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the probability parameters: n, m & h, and their voltage dependent rate constants (alpha and beta in [1/Sec]):


% n:
alpha_n = zeros(size(T), 'like', T);
beta_n = zeros(size(T), 'like', T);
n = zeros(size(T), 'like', T);
n(1) = n0;


% m:
alpha_m = zeros(size(T), 'like', T);
beta_m = zeros(size(T), 'like', T);
m = zeros(size(T), 'like', T);
m(1) = m0;


% h:
alpha_h = zeros(size(T), 'like', T);
beta_h = zeros(size(T), 'like', T);
h = zeros(size(T), 'like', T);
h(1) = h0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% voltage step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:(length(T)-1)
    
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
    % voltage:
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;

end


% a figure showing the existence of an action potential after a short pulse: 
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(4,1,1:3);
plot(T, V);
title('Voltage VS Time');
ylim([min(V)-10 max(V)+10]);                           
xlabel('Time [mS]');
ylabel({'Voltage', '[mV]'});
subplot(4,1,4);
plot(T, I_inj, Color = 'red');
title({'','Injected Current VS Time'});
xlabel('Time [mS]');
ylabel({'Current', 'I', '[nA]'});
ylim([min(I_inj)-10 max(I_inj)+10]);                           
sgtitle({'            an Action Potential due to a Short Pulse' , ''})
hold off;




%%   QUESTION 4   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% time and deltaT in [mS]:
deltaT = 0.01;
Time = 100;
T = 0:deltaT:Time;


% V_rest and voltage in [mV]:
V = zeros(size(T), 'like', T);
V(1) = V_rest;
dvdt = zeros(size(T), 'like', T);


% DEFINING REFRACTORY PERIOD (time after stimulus in mS)
REF_ABS_ms = 1;
REF_REL_ms = 3;
REF_ABS = REF_ABS_ms/deltaT;
REF_REL = REF_REL_ms/deltaT;


% I_inj (an external current injected to the nueron) in [nA]:         
I_inj = zeros(size(T), 'like', T);
CurrentInj = 15; 
CurrentDuration_ms = 0.5;
CurrentStartTime_ms = 10;
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = CurrentStartTime+CurrentDuration;
I_inj([CurrentStartTime:CurrentEndTime, CurrentStartTime+REF_ABS:CurrentEndTime+REF_ABS]) = CurrentInj;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n, m & h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the probability parameters: n, m & h, and their voltage dependent rate constants (alpha and beta in [1/Sec]):


% n:
alpha_n = zeros(size(T), 'like', T);
beta_n = zeros(size(T), 'like', T);
n = zeros(size(T), 'like', T);
n(1) = n0;


% m:
alpha_m = zeros(size(T), 'like', T);
beta_m = zeros(size(T), 'like', T);
m = zeros(size(T), 'like', T);
m(1) = m0;


% h:
alpha_h = zeros(size(T), 'like', T);
beta_h = zeros(size(T), 'like', T);
h = zeros(size(T), 'like', T);
h(1) = h0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% voltage step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:(length(T)-1)
    
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
    % voltage:
    dvdt(i) = (G_na*((m(i))^3)*h(i)*(E_na-V(i)) + G_k*((n(i))^4)*(E_k-V(i)) + G_m*(E_m-V(i)) + I_inj(i))/C_m;
    V(i+1) = V(i) + dvdt(i)*deltaT;

end


% a figure showing the existence of an action potential after a short pulse: 
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(4,1,1:3);
plot(T, V);
title('Voltage VS Time');
ylim([min(V)-10 max(V)+10]);                           
xlabel('Time [mS]');
ylabel({'Voltage', '[mV]'});
subplot(4,1,4);
plot(T, I_inj, Color = 'red');
title({'','Injected Current VS Time'});
xlabel('Time [mS]');
ylabel({'Current', 'I', '[nA]'});
ylim([min(I_inj)-10 max(I_inj)+10]);                           
sgtitle({'            an Action Potential due to a Short Pulse' , ''})
hold off;



