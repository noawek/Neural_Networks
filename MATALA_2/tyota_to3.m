
%%
figure('Units','normalized','position',[0 0 1 1]);
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [mS]')
ylabel('Current [nA]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:), "LineWidth", 2*neurons_num-2*i+1);
end
hold off
xlabel('Time [sec]')
ylabel('Rate [Hz]')
title('Rate vs. Time')
%ylim([-1.5 1.5])
legend(string(1:neurons_num))

%%

figure('Units','normalized','position',[0 0 1 1]);
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:), "LineWidth", 2*neurons_num-2*i+1);
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
ylim([-1 50])
legend(string(1:neurons_num))

%%
clear;
clc;



%%   QUESTION 3 (A)  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



neurons_num = 50;
W = zeros(neurons_num);
sigma1 = 30;
sigma2 = 60;
w_const = 0.2;


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-(w_const*exp(-((Di_j^2)/(sigma2^2))));
    end
end


% things related to time in [mS]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in Hz:
m = randperm(neurons_num);
n = 0.05:0.05:2.5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);


% initial random injected currents in [nA]:
I_inj = zeros(neurons_num,length(Time));
for j = 1:neurons_num
    I_inj(j,1:100) = randi(25,1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end


% plotting the results:
figure('Units','normalized','position',[0 0 1 1]);
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [mS]')
ylabel('Current [nA]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
hold off
xlabel('Time [mS]')
ylabel('Rate [Hz]')
title('Rate vs. Time')
ylim([-1 50])
legend(string(1:neurons_num))


%%
figure('Units','normalized','position',[0 0 1 1]);
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [sec]')
ylabel('Current [A]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:), "LineWidth", 2*neurons_num-2*i+1);
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
%ylim([-1.5 1.5])
legend(string(1:neurons_num))

%%

figure('Units','normalized','position',[0 0 1 1]);
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:), "LineWidth", 2*neurons_num-2*i+1);
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
ylim([-1 50])
legend(string(1:neurons_num))




%%   QUESTION 3 (B)  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons without inihibition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% no inihibition means that the ratio between the strength of inhibition
% and the strength of exitation (w_const), equals zero:
neurons_num = 50;
W = zeros(neurons_num);
sigma1 = 30;
sigma2 = 60;
w_const = 0;   


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-(w_const*exp(-((Di_j^2)/(sigma2^2))));
    end
end


% things related to time in [mS]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in Hz:
m = randperm(neurons_num);
n = 0.05:0.05:2.5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);


% initial random injected currents in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(1,1:100) = 0.7;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end


% plotting the results:
figure('Units','normalized','position',[0 0 1 1]);
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [sec]')
ylabel('Current [A]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
ylim([-1 50])
legend(string(1:neurons_num))





%%   QUESTION 1   %%   if i want to check per 1 degree each time




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading, ordering, and re-ordering the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



spiketrain_before = 0;
APs = 0;
APs1 = 0;
APs2 = 0;
AP_perm1 = 0;
AP_perm2 = 0;
S = dir(fullfile('*.mat'));


for k = 1:numel(S)
    

    F = fullfile(S(k).name);
    load(F);
    head_direction_deg = ((headDirection(:))')*180/pi;
    
    for i = 1:360
        deg(i,:) = [(i-1), (i)];
    end
    
    
    if length(spiketrain_before) == length(spiketrain)
         
        for n = 1:length(spiketrain1)
             
            if spiketrain1(n) == 1 && APs1(1) == 0 
                 APs1 = n;
             elseif spiketrain1(n) == 1 && APs1(1) ~= 0
                 APs1(length(APs1)+1) = n - APs1_before;
             end
             
             if spiketrain1(n) == 1
                 APs1_before = n;
             end

         end
         
         for l = 1:length(spiketrain2)
             
             if spiketrain2(l) == 1 && APs2(1) == 0 
                 APs2 = l;
             elseif spiketrain2(l) == 1 && APs2(1) ~= 0
                 APs2(length(APs2)+1) = l - APs2_before;
             end
             
             if spiketrain2(l) == 1
                 APs2_before = l;
             end

         end
         
         AP_perm1 = APs1(randperm(length(APs1)));
         spiketrain_fake = zeros(1,length(spiketrain1));
         num = 0;
         
         for a = 1:length(AP_perm1)
             num = num+AP_perm1(a);
             spiketrain_fake(num) = 1;
         end

         spiketrain_fakes{k} = spiketrain_fake;
         
         AP_perm2 = APs2(randperm(length(APs2)));
         spiketrain_fake = zeros(1,length(spiketrain2));
         num = 0;

         for c = 1:length(AP_perm2)
             num = num+AP_perm2(c);
             spiketrain_fake(num) = 1;
         end

         spiketrain_fakes{k+1} = spiketrain_fake;
   
        for j = 1:length(deg)
             places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
             countperdeg = length(places); 
             timeperdeg = countperdeg*(1/sampleRate);
             num_spikes1 = sum(spiketrain1(places));
             num_spikes2 = sum(spiketrain2(places));
             rateperdeg(k,j) = num_spikes1/timeperdeg;
             rateperdeg(k+1,j) = num_spikes2/timeperdeg;
             num_spikes1_fake = sum(spiketrain_fakes{k}(places));
             num_spikes2_fake = sum(spiketrain_fakes{k+1}(places));
             rateperdeg_fake(k,j) = num_spikes1_fake/timeperdeg;
             rateperdeg_fake(k+1,j) = num_spikes2_fake/timeperdeg;
         end
         
    else
   
    APs = 0;

    for m = 1:length(spiketrain)
        
        if spiketrain(m) == 1 && APs(1) == 0 
            APs = m;
        elseif spiketrain(m) == 1 && APs(1) ~= 0
            APs(length(APs)+1) = m - APs_before;
        end

        if spiketrain(m) == 1
                 APs_before = m;
        end

    end

    AP_perm = 0;
    AP_perm = APs(randperm(length(APs)));
    spiketrain_fake = zeros(1,length(spiketrain));
    num = 0;

    for b = 1:length(AP_perm)
        num = num+AP_perm(b);
        spiketrain_fake(num) = 1;
    end

    spiketrain_fakes{k} = spiketrain_fake;

    for j = 1:length(deg)
        places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
        countperdeg = length(places); 
        timeperdeg = countperdeg*(1/sampleRate);
        num_spikes = sum(spiketrain(places));
        num_spikes_fake = sum(spiketrain_fake(places));
        rateperdeg(k,j) = num_spikes/timeperdeg;
        rateperdeg_fake(k,j) = num_spikes_fake/timeperdeg;
    end

    end

    C{k} = load(F);
    spiketrain_before = spiketrain;
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:size(rateperdeg,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(1,2,1);
    plot([deg(:,1);360],[rateperdeg(i,:),rateperdeg(i,1)]);
    xlim([0,360]);
    ylim([0, 17]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title('Line Graph');
    subplot(1,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg(i,:),rateperdeg(i,1)]);
    rlim([0, 17]);
    title('Polar Plot');
    legend('Firing Rate');
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rate [Hz] as a function of the Head Direction [degrees]'});
    hold off;
end


clear S F AP_perm1 AP_perm2 APs APs1 APs2 APs2_before APs1_before APs_before boxSize countperdeg head_direction_deg headDirection sampleRate spiketrain1 spiketrain2 spiketrain_fake spiketrain AP_perm ans j b m k c l a n i num num_spikes2_fake num_spikes1_fake num_spikes num_spikes_fake num_spikes2 num_spikes1 spiketrain_before timeperdeg posx posy post theta phase places;














%%   QUESTION 3   %%




% neurons:

neurons_num = 50;
W = zeros(neurons_num);
sigma1 = 30;
sigma2 = 60;
w_const = 0.2;

for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-(w_const*exp(-((Di_j^2)/(sigma2^2))));
    end
end


% things related to time in [mS]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;


% rate in Hz:
m = randperm(neurons_num);
n = 0.05:0.05:2.5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);
%%
% current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(1,1:100) = 30;
I_inj(3,1:100) = 22;



for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end


figure('Units','normalized','position',[0 0 1 1]);
subplot(2,1,1)
hold on
plot(Time,I_inj)
xlabel('Time [sec]')
ylabel('Current [A]')
title('Current vs. Time')
subplot(2,1,2)
hold on
for i = 1:neurons_num 
    plot(Time,rate(i,:))
end
hold off
xlabel('Time [sec]')
ylabel('Rate')
title('Rate vs. Time')
ylim([-1 40])
legend(string(1:neurons_num))




%%
for T = 1:1000



end



%%
x = -2:2;
xp = subplus(x);


plot(x,xp)
ylim([-0.5 2.5])

%%
figure
plot(Time,p);
ylim([min(p) max(p)]);

p = rate(1,:);

%%   QUESTION 3 (A)  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
neurons_num = 50;
W = zeros(neurons_num);
sigma1 = 30;
sigma2 = 60;
w_const = 0.2;


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-(w_const*exp(-((Di_j^2)/(sigma2^2))));
    end
end


% things related to time in [sec]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;


% no injected currents in [nA]:
I_inj = zeros(neurons_num,length(Time));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.1:0.1:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end


% plotting the results until 50 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(3,1,1);
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title('until 50 Hz');
ylim([-1 50]);
subplot(3,1,2);
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title('until 100 Hz');
ylim([-1 100]);
subplot(3,1,3);
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title('until 150 Hz');
ylim([-1 150]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
sgtitle({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
hold off;


