



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear;
clc;





%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading, ordering, and re-ordering the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



spiketrain_before = 0;
APs = 0;
APs1 = 0;
APs2 = 0;
AP_perm1 = 0;
AP_perm2 = 0;
AP_permute1 = 0;
AP_permute2 = 0;
S = dir(fullfile('*.mat'));



for k = 1:numel(S)

    
    % loading the files and the relevant data from them:
    F = fullfile(S(k).name);
    load(F);
    head_direction_deg = ((headDirection(:))')*180/pi;
    

    % setting the bins of degrees for the analyses (each 10 are grouped together):
    for i = 1:36
        deg(i,:) = [((i-1)*10), (i*10)];
    end
    
    
    % extracting some important variables from the files (using "if" to discriminate the first 4 files, each includes data only on one
    % neuron recorded, from the last file that includes data recorded from 2 neurons simultaniously):  

   
    % extracting spike trains of the two last neurons (in the last file), so we can calculate rates later on:
    if length(spiketrain_before) == length(spiketrain)        
        
        % neuron 5:
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
        
        % neuron 6:
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
         

         % permuting the spike trains data as a preperation for Question 2: 

         % first permutation - neuron 5:
         AP_perm1 = APs1(randperm(length(APs1)));
         spiketrain_fake = zeros(1,length(spiketrain1));
         num = 0;     
         for a = 1:length(AP_perm1)
             num = num+AP_perm1(a);
             spiketrain_fake(num) = 1;
         end
         spiketrain_fakes{k} = spiketrain_fake;
         
         % first permutation - neuron 6:
         AP_perm2 = APs2(randperm(length(APs2)));
         spiketrain_fake = zeros(1,length(spiketrain2));
         num = 0;
         for c = 1:length(AP_perm2)
             num = num+AP_perm2(c);
             spiketrain_fake(num) = 1;
         end
         spiketrain_fakes{k+1} = spiketrain_fake;
         
         % second permutation - neuron 5:
         AP_permute1 = APs1(randperm(length(APs1)));
         spiketrain_fake = zeros(1,length(spiketrain1));
         num = 0;         
         for a = 1:length(AP_permute1)
             num = num+AP_permute1(a);
             spiketrain_fake(num) = 1;
         end
         spiketrain_fakes2{k} = spiketrain_fake;
         
         % second permutation - neuron 6:
         AP_permute2 = APs2(randperm(length(APs2)));
         spiketrain_fake = zeros(1,length(spiketrain2));
         num = 0;
         for c = 1:length(AP_permute2)
             num = num+AP_permute2(c);
             spiketrain_fake(num) = 1;
         end
         spiketrain_fakes2{k+1} = spiketrain_fake;


         % calculation of rates (2 fakes and 1 original) for neurons 5 + 6:
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
             num_spikes1_fake2 = sum(spiketrain_fakes2{k}(places));
             num_spikes2_fake2 = sum(spiketrain_fakes2{k+1}(places));
             rateperdeg_fake2(k,j) = num_spikes1_fake2/timeperdeg;
             rateperdeg_fake2(k+1,j) = num_spikes2_fake2/timeperdeg;
         end
    
         
    % extracting the same important variables and preforming relevant calculations, 
    % but this time for the rest of the neurons (1-4):
    else


    % zeroing the APs variable each the the loop runs (so data won't mix):
    APs = 0;
    
    % extracting spike trains of those neurons (from the first 4 files), so we can calculate rates later on:
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


    % permuting the spike trains data as a preperation for Question 2:

    % first permutation:
    AP_perm = 0;
    AP_perm = APs(randperm(length(APs)));
    spiketrain_fake = zeros(1,length(spiketrain));
    num = 0;
    for b = 1:length(AP_perm)
        num = num+AP_perm(b);
        spiketrain_fake(num) = 1;
    end
    spiketrain_fakes{k} = spiketrain_fake;
    
    % second permutation:
    AP_permute = 0;
    AP_permute = APs(randperm(length(APs)));
    spiketrain_fake = zeros(1,length(spiketrain));
    num = 0;
    for b = 1:length(AP_permute)
        num = num+AP_permute(b);
        spiketrain_fake(num) = 1;
    end
    spiketrain_fakes2{k} = spiketrain_fake;


    % calculation of rates (2 fakes and 1 original) for neurons all 4 neurons:
    for j = 1:length(deg)
        places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
        countperdeg = length(places); 
        timeperdeg = countperdeg*(1/sampleRate);
        num_spikes = sum(spiketrain(places));
        num_spikes_fake = sum(spiketrain_fakes{k}(places));
        num_spikes_fake2 = sum(spiketrain_fakes2{k}(places));
        rateperdeg(k,j) = num_spikes/timeperdeg;
        rateperdeg_fake(k,j) = num_spikes_fake/timeperdeg;
        rateperdeg_fake2(k,j) = num_spikes_fake2/timeperdeg;
    end


    end


    % just so we'll have the data in case we need it...
    C{k} = load(F);

   
    % a way to track whether we got to the last file (which contains data
    % about 2 neurons) or not (so if we did get there, it will activate the
    % "if" loop to act in a certain way):
    spiketrain_before = spiketrain;
    

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures of all 6 neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:size(rateperdeg,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(1,2,1);
    plot([deg(:,1);360],[rateperdeg(i,:),rateperdeg(i,1)]);
    xlim([0,360]);
    % 14 Hz is the maximum firing rate of all neurons, and setting this limit is meant for us to be able to compare different neurons:
    ylim([0, 14]);    
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title('Line Graph');
    subplot(1,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg(i,:),rateperdeg(i,1)]);
    % 14 Hz again for the same reason:
    rlim([0, 14]);    
    title('Polar Plot');
    legend('Firing Rate');
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rate [Hz] as a function of the Head Direction [degrees]'});
    hold off;
end


% clearing some irrelevant variables so the workspace would be comfetrable to work with:
clear S F AP_perm1 AP_perm2 APs APs1 APs2 APs2_before APs1_before APs_before boxSize countperdeg head_direction_deg headDirection sampleRate spiketrain1 spiketrain2 spiketrain_fake spiketrain AP_perm ans j b m k c l a n i num num_spikes2_fake num_spikes1_fake num_spikes num_spikes_fake num_spikes2 num_spikes1 spiketrain_before timeperdeg posx posy post theta phase places;





%%   QUESTION 2   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% creating the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% we already managed to have a vector of the margines between action potentials.
% let's play with it now by permuting it:

for i = 1:size(rateperdeg_fake,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(1,2,1);
    plot([deg(:,1);360],[rateperdeg_fake(i,:),rateperdeg_fake(i,1)]);
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title('Line Graph');
    subplot(1,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg_fake(i,:),rateperdeg_fake(i,1)]);
    rlim([0, 14]);
    title('Polar Plot');
    legend('Firing Rate');
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rate [Hz] as a function of the Head Direction [degrees]', '\fontsize{10} (the inter-spike interval of the spike train was randomly permuted)'});
    hold off;
end


clear i


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comparing the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:size(rateperdeg_fake,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(2,2,1);
    plot([deg(:,1);360],[rateperdeg(i,:),rateperdeg(i,1)]);
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title({'Line Graph', 'Original Data'});
    subplot(2,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg(i,:),rateperdeg(i,1)]);
    rlim([0, 14]);
    title({'Polar Plot', 'Original Data'});
    legend('Firing Rate');
    subplot(2,2,3);
    plot([deg(:,1);360],[rateperdeg_fake(i,:),rateperdeg_fake(i,1)], Color = 'r');
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title({'Line Graph', 'Randomly Permuted Data'});
    subplot(2,2,4);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg_fake(i,:),rateperdeg_fake(i,1)], Color = 'r');
    rlim([0, 14]);
    title({'Polar Plot', 'Randomly Permuted Data'});
    legend('Firing Rate');
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rates [Hz] as a function of Head Direction [degrees]:', 'comparison between the original data from the neuron, and the test data (where the inter-spike interval of the spike train was randomly permuted)', '\fontsize{10} to check if a neuron can really be classified as a head direction cell'});
    hold off;
end


clear i



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comparing the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:size(rateperdeg_fake,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(3,2,1);
    plot([deg(:,1);360],[rateperdeg(i,:),rateperdeg(i,1)]);
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title({'Line Graph', 'Original Data'});
    subplot(3,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg(i,:),rateperdeg(i,1)]);
    rlim([0, 14]);
    title({'Polar Plot', 'Original Data'});
    legend('Firing Rate');
    subplot(3,2,3);
    plot([deg(:,1);360],[rateperdeg_fake(i,:),rateperdeg_fake(i,1)], Color = 'r');
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title({'Line Graph', 'Randomly Permuted Data'});
    subplot(3,2,4);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg_fake(i,:),rateperdeg_fake(i,1)], Color = 'r');
    rlim([0, 14]);
    title({'Polar Plot', 'Randomly Permuted Data'});
    legend('Firing Rate');
    subplot(3,2,5);
    plot([deg(:,1);360],[rateperdeg_fake2(i,:),rateperdeg_fake2(i,1)], Color = 'g');
    xlim([0,360]);
    ylim([0, 14]);
    ylabel('Firing Rate [Hz]');
    xlabel('Head Direction in Degrees [\circ]');
    title({'Line Graph', 'Randomly Permuted Data'});
    subplot(3,2,6);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg_fake2(i,:),rateperdeg_fake2(i,1)], Color = 'g');
    rlim([0, 14]);
    title({'Polar Plot', 'Randomly Permuted Data'});
    legend('Firing Rate');
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rates [Hz] as a function of Head Direction [degrees]:', 'comparison between the original data from the neuron, and the test data (where the inter-spike interval of the spike train was randomly permuted)', '\fontsize{10} to check if a neuron can really be classified as a head direction cell'});
    hold off;
end


clear i




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







%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%



clear;
clc;




%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


spiketrain_before = 0;
S = dir(fullfile('*.mat'));
for k = 1:numel(S)
    F = fullfile(S(k).name);
    load(F);
    head_direction_deg = ((headDirection(:))')*180/pi;
    for i = 1:36
        deg(i,:) = [((i-1)*10), (i*10)];
    end
    if length(spiketrain_before) == length(spiketrain)
         for j = 1:length(deg)
             places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
             countperdeg = length(places); 
             timeperdeg = countperdeg*0.001;
             num_spikes1 = sum(spiketrain1(places));
             num_spikes2 = sum(spiketrain2(places));
             rateperdeg(k,j) = num_spikes1/timeperdeg;
             rateperdeg(k+1,j) = num_spikes2/timeperdeg;
         end
    else
    for j = 1:length(deg)
        places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
        countperdeg = length(places); 
        timeperdeg = countperdeg*0.001;
        num_spikes = sum(spiketrain(places));
        rateperdeg(k,j) = num_spikes/timeperdeg;
    end
    end
    C{k} = load(F);
    spiketrain_before = spiketrain;
end
%%


for i = 1:size(rateperdeg,1)

    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    subplot(1,2,1);
    plot([deg(:,1);360],[rateperdeg(i,:),rateperdeg(i,1)]);
    xlim([0,360]);
    ylim([0, 14]);
    subplot(1,2,2);
    polarplot(deg2rad([deg(:,1);360]), [rateperdeg(i,:),rateperdeg(i,1)]);
    rlim([0, 14]);
    hold off;


end


%% neuron 1:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot([deg(:,1);360],[rateperdeg(6,:),rateperdeg(6,1)]);
hold off;

%%
[theta, rho] = cart2pol(deg(:,2),rateperdeg(6,:));
figure('Units','normalized','position',[0 0 1 1]);

polarplot(deg2rad([deg(:,1);360]), [rateperdeg(6,:),rateperdeg(6,1)]);

%%

S = dir(fullfile('*.mat'));
for k = 1:numel(S)
    F = fullfile(S(k).name);
    load(F);
    head_direction_deg = ((headDirection(:))')*180/pi;
    for i = 1:72
        deg(i,:) = [((i-1)*10), (i*10)];
    end
    if length(spiketrain1) == length(spiketrain2)
         for j = 1:length(deg)
             places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
             countperdeg(k,j,:) = length(places); 
             timeperdeg(k,j,:) = countperdeg(k,j,:)*0.02;
             num_spikes(k,j,:) = sum(spiketrain1(places));
             num_spikes(k+1,j,:) = sum(spiketrain2(places));
             rateperdeg(k,j,:) = num_spikes(j)/timeperdeg(j);
             plot(deg(:,2),rateperdeg(j,:));
         end
    else
    for j = 1:length(deg)
        places = find(head_direction_deg < deg(j,2) & head_direction_deg > deg(j,1));
        countperdeg(j,:) = length(places); 
        timeperdeg(j,:) = countperdeg(j)*0.02;
        num_spikes(j,:) = sum(spiketrain(places));
        rateperdeg(j,:) = num_spikes(j)/timeperdeg(j);
        plot(deg(:,2),rateperdeg(j,:));


    end
    end
    C{k} = load(F);
end


%%
S = dir(fullfile('*.mat'));
for k = 1:numel(S)
    F = fullfile(S(k).name);
    load(F);
    head_direction{k} = (headDirection(:))';
    C{k} = load(F);
end

%%
if AP{1} == AP{5}
    disp 'yay'
else
    disp 'nay'
end
%%
unequal = 0;
for i = 1:length(AP{5})

    if AP{5}(i) ~= AP{1}(i)

        unequal = [unequal,i];

    end

end