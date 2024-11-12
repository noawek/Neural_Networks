%%   QUESTION 4   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settimg some parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% time and deltaT in [mS]:
deltaT = 0.01;
Time = 50;
T = 0:deltaT:Time;


% V_rest and voltage in [mV]:
V = zeros(size(T), 'like', T);
V(1) = V_rest;
dvdt = zeros(size(T), 'like', T);


% DEFINING REFRACTORY PERIOD (time after stimulus in mS)
REF_ABS_ms = 1;
REF_REL_ms = 10;
REF_ABS = REF_ABS_ms/deltaT;
REF_REL = REF_REL_ms/deltaT;


% I_inj (an external current injected to the nueron) in [nA]:         
I_inj = zeros(size(T), 'like', T);
CurrentInj1 = 15;
CurrentInj2 = 70; 
CurrentDuration_ms = 0.5;
CurrentStartTime_ms = 10;
CurrentDuration = CurrentDuration_ms/deltaT;
CurrentStartTime = (CurrentStartTime_ms/deltaT)+1;
CurrentEndTime = CurrentStartTime+CurrentDuration;
I_inj(CurrentStartTime:CurrentEndTime) = CurrentInj1;
I_inj(CurrentStartTime+REF_REL:CurrentEndTime+REF_REL) = CurrentInj2;




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


