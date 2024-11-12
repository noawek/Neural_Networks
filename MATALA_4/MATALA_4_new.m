



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
clc;





%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
neurons_num = 200;
W = zeros(neurons_num);
sigma1 = 10;


% the W's matrix (setting the synapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% choosing one neuron for all demonstrations:
the_neuron = 100;


% 3 neurons for the upcoming demonstration (including "the_neuron"):
some_neurons  = [(the_neuron-50),the_neuron,(the_neuron+50)];


% plotting the synaptic weights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(2,3,1);
hold on;
plot((0:neurons_num),[W(200,some_neurons(1));W(:,some_neurons(1))]);
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm considering the neuron numbers'});
subplot(2,3,2);
hold on;
plot((0:neurons_num),[W(200,some_neurons(2));W(:,some_neurons(2))], Color='g');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm considering the neuron numbers'});
subplot(2,3,3);
hold on;
plot((0:neurons_num),[W(200,some_neurons(3));W(:,some_neurons(3))], Color='m');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm considering the neuron numbers'});
subplot(2,3,4);
hold on;
plot((0:neurons_num)-some_neurons(1),[W(200,some_neurons(1));W(:,some_neurons(1))]);
plot((0:neurons_num)-some_neurons(1),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('the Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,5);
hold on;
plot((0:neurons_num)-some_neurons(2),[W(200,some_neurons(2));W(:,some_neurons(2))], Color='g');
plot((0:neurons_num)-some_neurons(2),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,6);
hold on;
plot((0:neurons_num)-some_neurons(3),[W(200,some_neurons(3));W(:,some_neurons(3))], Color='m');
plot((0:neurons_num)-some_neurons(3),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm by the distance of neurons from this neuron'});
sgtitle({'\fontsize{14} \bf Synaptic Weights (W) vs. the Distance from a Neuron', '\fontsize{12} \rm 3 examples of chosen neurons', 'from 200 neurons affecting each other in a ring structure',''});
hold off;


% things related to time in [sec]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;


% initial injected current to the chosen neuron in [nA]:
I_intensity = 50;
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = I_intensity;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end


% normalizing the rates: 
rate(:,1) = rate(:,1)/max(rate(:,1));
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('center_neurons.avi');
open(vidfile);
for i = 1:(length(Time))
    % a movie of the model:
    plot((0:neurons_num),[R(200);R]);
    xlabel('Neuron Number');
    ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
    title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 200 neurons affecting each other in a ring structure', 'starting with different random firing rates', ['\fontsize{12} \rm while neuron number ' num2str(the_neuron) ' gets a current of a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning'], ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec']});
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    % normalizing:
    R = R/max(R);      
    rate(:,i+1) = R;
end
close(vidfile);
hold off;


% fixing the rate:
rate = rate(:,1:(size(rate,2)-1));


% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
title({'\fontsize{14} Rate vs. Time', ['\fontsize{12} \rm of neuron number ' num2str(the_neuron)], 'who is the only one (from 200 neurons in a ring model, starting with random firing rates)', ['that got a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning']});
hold off;





%%   QUESTION 2 (A) - creating a ring with a diversion to the left  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_l = zeros(neurons_num);


% the W's matrix (setting the synapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_l,1)
    for j = 1:size(W_l,2)
        Di_j = min(abs((j-2)-i), neurons_num - abs((j-2)-i));
        W_l(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% plotting the synaptic weights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(2,3,1);
hold on;
plot((0:neurons_num),[W_l(200,some_neurons(1));W_l(:,some_neurons(1))]);
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm considering the neuron numbers'});
subplot(2,3,2);
hold on;
plot((0:neurons_num),[W_l(200,some_neurons(2));W_l(:,some_neurons(2))], Color='g');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm considering the neuron numbers'});
subplot(2,3,3);
hold on;
plot((0:neurons_num),[W_l(200,some_neurons(3));W_l(:,some_neurons(3))], Color='m');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm considering the neuron numbers'});
subplot(2,3,4);
hold on;
plot((0:neurons_num)-some_neurons(1),[W_l(200,some_neurons(1));W_l(:,some_neurons(1))]);
plot((0:neurons_num)-some_neurons(1),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('the Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,5);
hold on;
plot((0:neurons_num)-some_neurons(2),[W_l(200,some_neurons(2));W_l(:,some_neurons(2))], Color='g');
plot((0:neurons_num)-some_neurons(2),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,6);
hold on;
plot((0:neurons_num)-some_neurons(3),[W_l(200,some_neurons(3));W_l(:,some_neurons(3))], Color='m');
plot((0:neurons_num)-some_neurons(3),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm by the distance of neurons from this neuron'});
sgtitle({'\fontsize{12} Counterclockwise cells', '\fontsize{14} \bf Synaptic Weights (W) vs. the Distance from a Neuron', '\fontsize{12} \rm 3 examples of chosen neurons', 'from 200 neurons affecting each other in a ring structure',''});
hold off;


% initial injected current to the chosen neuron in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = I_intensity;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end


% normalizing the rates: 
rate(:,1) = rate(:,1)/max(rate(:,1));
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('left_neurons.avi');
open(vidfile);
for i = 1:(length(Time))
    % a movie of the model:
    plot((0:neurons_num),[R(200);R]);
    xlabel('Neuron Number');
    ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
    title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm (Counterclockwise cells)', 'of 200 neurons affecting each other in a ring structure', 'starting with different random firing rates', ['\fontsize{12} \rm while neuron number ' num2str(the_neuron) ' gets a current of a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning'], ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec']});
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W_l*R));
    % normalizing:
    R = R/max(R);
    rate(:,i+1) = R;
end
close(vidfile);
hold off;


% fixing the rate:
rate = rate(:,1:(size(rate,2)-1));


% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
title({'\fontsize{14} Rate vs. Time', ['\fontsize{12} \rm of neuron number ' num2str(the_neuron) ' (from the Counterclockwise cells)'], 'who is the only one (from 200 neurons in a ring model, starting with random firing rates)', ['that got a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning']});
hold off;





%%   QUESTION 2 (B) - creating a ring with a diversion to the right  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_r = zeros(neurons_num);


% the W's matrix (setting the synapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_r,1)
    for j = 1:size(W_r,2)
        Di_j = min(abs((j+2)-i), neurons_num - abs((j+2)-i));
        W_r(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% plotting the synaptic weights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
subplot(2,3,1);
hold on;
plot((0:neurons_num),[W_r(200,some_neurons(1));W_r(:,some_neurons(1))]);
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm considering the neuron numbers'});
subplot(2,3,2);
hold on;
plot((0:neurons_num),[W_r(200,some_neurons(2));W_r(:,some_neurons(2))], Color='g');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm considering the neuron numbers'});
subplot(2,3,3);
hold on;
plot((0:neurons_num),[W_r(200,some_neurons(3));W_r(:,some_neurons(3))], Color='m');
plot(0:neurons_num,zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Neuron Number');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm considering the neuron numbers'});
subplot(2,3,4);
hold on;
plot((0:neurons_num)-some_neurons(1),[W_r(200,some_neurons(1));W_r(:,some_neurons(1))]);
plot((0:neurons_num)-some_neurons(1),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('the Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(1))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,5);
hold on;
plot((0:neurons_num)-some_neurons(2),[W_r(200,some_neurons(2));W_r(:,some_neurons(2))], Color='g');
plot((0:neurons_num)-some_neurons(2),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(2))]; '\rm by the distance of neurons from this neuron'});
subplot(2,3,6);
hold on;
plot((0:neurons_num)-some_neurons(3),[W_r(200,some_neurons(3));W_r(:,some_neurons(3))], Color='m');
plot((0:neurons_num)-some_neurons(3),zeros(neurons_num+1,1),LineStyle="--", Color='r');
xlabel('Distance from the neuron');
ylabel('Synaptic Weight [W]');
title({['Synaptic Weights [W] of neuron number ' num2str(some_neurons(3))]; '\rm by the distance of neurons from this neuron'});
sgtitle({'\fontsize{12} Clockwise cells', '\fontsize{14} \bf Synaptic Weights (W) vs. the Distance from a Neuron', '\fontsize{12} \rm 3 examples of chosen neurons', 'from 200 neurons affecting each other in a ring structure',''});
hold off;


% initial injected current to the chosen neuron in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = I_intensity;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end


% normalizing the rates: 
rate(:,1) = rate(:,1)/max(rate(:,1));
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('right_neurons.avi');
open(vidfile);
for i = 1:(length(Time))
    % a movie of the model:
    plot((0:neurons_num),[R(200);R]);
    xlabel('Neuron Number');
    ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
    title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm (Clockwise cells)', 'of 200 neurons affecting each other in a ring structure', 'starting with different random firing rates', ['\fontsize{12} \rm while neuron number ' num2str(the_neuron) ' gets a current of a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning'], ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec']});
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W_r*R));
    % normalizing:
    R = R/max(R);
    rate(:,i+1) = R;
end
close(vidfile);
hold off;


% fixing the rate:
rate = rate(:,1:(size(rate,2)-1));


% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in that point in time]'});
title({'\fontsize{14} Rate vs. Time', ['\fontsize{12} \rm of neuron number ' num2str(the_neuron) ' (from the Clockwise cells)'], 'who is the only one (from 200 neurons in a ring model, starting with random firing rates)', ['that got a current of ' num2str(I_inj(the_neuron, 1)) 'nA  for ' num2str(10*dt) 'sec at the beginning']});
hold off;





%%   QUESTION 3   %% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% injected current in [nA]:
I_inj_l = zeros(neurons_num,length(Time));
I_inj_l(:, 200:400) = I_intensity;
I_inj_r = zeros(neurons_num,length(Time));
I_inj_r(:, 600:800) = I_intensity;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;


% centered neurons:
rate_c = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_c(k,1) = n(m(k));
end
% normalizing the rates: 
rate_c(:,1) = rate_c(:,1)/max(rate_c(:,1));
R_c = rate_c(:,1);


% counterclockwise neurons:
rate_l = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_l(k,1) = n(m(k));
end
% normalizing the rates: 
rate_l(:,1) = rate_l(:,1)/max(rate_l(:,1));
R_l = rate_l(:,1);


% clockwise neurons:
rate_r = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_r(k,1) = n(m(k));
end
% normalizing the rates: 
rate_r(:,1) = rate_r(:,1)/max(rate_r(:,1));
R_r = rate_r(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('all_neurons.avi');
open(vidfile);
for i = 1:(length(Time))
    % a movie of the model:
    subplot(3,1,1);
    plot((0:neurons_num),[rate_c(200,i);rate_c(:,i)]);
    xlabel('Neuron Number');
    title({'\fontsize{12} Pure cells (main ring):', '\fontsize{10} \rm 200 neurons affecting each other (and affected by neurons from 2 other rings) in a ring structure'});
    subplot(3,1,2);
    plot((0:neurons_num),[rate_l(200,i);rate_l(:,i)], Color='r');
    ylim([min(rate_l(:,i)) 1]);
    xlabel('Neuron Number');
    ylabel({'Rate (norm)'; '[Rate in the time mentioned [Hz]/the maximum value from all neurons Rates in this ring in that point of time]'});
    title({'', '\fontsize{12} Counterclockwise cells:', '\fontsize{10} \rm 200 neurons affecting each other (and affected by neurons from 2 other rings) in a ring structure'});
    subplot(3,1,3);
    plot((0:neurons_num),[rate_r(200,i);rate_r(:,i)], Color='g');
    ylim([min(rate_r(:,i)) 1]);
    xlabel('Neuron Number');
    title({'', '\fontsize{12} Clockwise cells:', '\fontsize{10} \rm 200 neurons affecting each other (and affected by neurons from 2 other rings) in a ring structure'});
    if I_inj_l(1,i) == I_intensity
        sgtitle({'\fontsize{14} \bf Rate vs. Time', '\fontsize{12} \rm of 600 neurons affecting each other in 3 different ring structures, starting with random firing rates', ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec'], '\color{red} Head rotation to the left (current injected to the ring of the counterclockwise cells)'});
    elseif I_inj_r(1,i) == I_intensity
        sgtitle({'\fontsize{14} \bf Rate vs. Time', '\fontsize{12} \rm of 600 neurons affecting each other in 3 different ring structures, starting with random firing rates', ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec'], '\color{green} Head rotation to the right (current injected to the ring of the clockwise cells)'});
    else
        sgtitle({'\fontsize{14} \bf Rate vs. Time', '\fontsize{12} \rm of 600 neurons affecting each other in 3 different ring structures, starting with random firing rates', ['\bf Time = ' num2str(floor((i-1)*dt)) ' sec']});
    end
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R_c = rate_c(:,i) + dt*(1/Tau)*(-rate_c(:,i) + subplus(W*rate_c(:,i) + W*rate_l(:,i) + W*rate_r(:,i)));
    R_l = rate_l(:,i) + dt*(1/Tau)*(-rate_l(:,i) + subplus(W_l*rate_l(:,i) + W*rate_c(:,i) + W*rate_r(:,i)) + I_inj_l(:,i));
    R_r = rate_r(:,i) + dt*(1/Tau)*(-rate_r(:,i) + subplus(W_r*rate_r(:,i) + W*rate_c(:,i) + W*rate_l(:,i)) + I_inj_r(:,i));
    % normalizing: 
    R_c = R_c/max(R_c);
    rate_c(:,i+1) = R_c;
    R_l = R_l/max(R_l);
    rate_l(:,i+1) = R_l;
    R_r = R_r/max(R_r);
    rate_r(:,i+1) = R_r;
end
close(vidfile);
hold off;


% fixing the rate:
rate_c = rate_c(:,1:(size(rate_c,2)-1));
rate_l = rate_l(:,1:(size(rate_l,2)-1));
rate_r = rate_r(:,1:(size(rate_r,2)-1));




