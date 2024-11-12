


%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%



clear;
clc;




%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = dir(fullfile('*.mat'));
for k = 1:numel(S)
    F = fullfile(S(k).name);
    head_direction(k,:) = headDirection;
    C{k} = load(F);
end


%%


S = dir(fullfile('*.mat'));
for k = 1:numel(S)
    F = fullfile(S(k).name);
    head_direction(k,:) = headDirection;
    C{k} = load(F);
    for i = 1:length(c{k=='headDirection'})
        (headDirection_num2str(k)) = 


    end


end

for i = 







%%   QUESTION 3   %%




% neurons:

neurons_num = 50;
W = zeros(neurons_num);
sigma1 = 30;
sigma2 = 60;
w_const = 0.2;

for i = 1:size(W,1)
    for j = 1:size(W,2)
        W(i,j) = exp(-((abs(j-i)^2)/(sigma1^2)))-(w_const*exp(-((abs(j-i)^2)/(sigma2^2))));
    end
end


% things related to time in [mS]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;


% rate in Hz:
r = zeros(neurons_num,length(Time));
R = zeros(neurons_num,1);

% current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(1,1:100) = 0.7;


for i = 1:length(Time)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    r(:,i) = R;
end


figure
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [sec]')
ylabel('Current [A]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(1:length(Time),r(i,:))
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
ylim([-1.5 1.5])
legend(string(1:neurons_num))




%%
for T = 1:1000



end



%%
x = -2:2;
xp = subplus(x);


plot(x,xp)
ylim([-0.5 2.5])










