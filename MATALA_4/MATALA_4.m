



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
clc;





%%   QUESTION 1   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
neurons_num = 200;
W = zeros(neurons_num);
sigma1 = 10;


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W,1)
    for j = 1:size(W,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% choosing one neuron for all demonstrations:
the_neuron = 100;


% plotting the sinaptic wheights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot((1:neurons_num),W(:,the_neuron));
plot(1:neurons_num,zeros(neurons_num,1),LineStyle="--",Color='r');
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;


% things related to time in [sec]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;


% initial injected current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('center_neurons.avi');
open(vidfile);
for i = 1:(length(Time)-1)
    % a movie of the model:
    plot(R/max(R));
    xlabel('neurons');
    ylabel('Rate');
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end
close(vidfile);
hold off;


% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;





%%   QUESTION 2 (A) - creating a ring with a diversion to the left  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_l = zeros(neurons_num);


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_l,1)
    for j = 1:size(W_l,2)
        Di_j = min(abs((j-2)-i), neurons_num - abs((j-2)-i));
        W_l(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% plotting the sinaptic wheights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot((1:neurons_num),W_l(:,the_neuron));
plot(1:neurons_num,zeros(neurons_num,1),LineStyle="--",Color='r');
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;


% initial injected current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('left_neurons.avi');
open(vidfile);
for i = 1:(length(Time)-1)
    % a movie of the model:
    plot(R/max(R));
    xlabel('neurons');
    ylabel('Rate');
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W_l*R));
    rate(:,i+1) = R;
end
close(vidfile);
hold off;

% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;





%%   QUESTION 2 (B) - creating a ring with a diversion to the right  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_r = zeros(neurons_num);


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_r,1)
    for j = 1:size(W_r,2)
        Di_j = min(abs((j+2)-i), neurons_num - abs((j+2)-i));
        W_r(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end


% plotting the sinaptic wheights (W) vs. the distance from a neuron:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot((1:neurons_num),W_r(:,the_neuron));
plot(1:neurons_num,zeros(neurons_num,1),LineStyle="--",Color='r');
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;


% initial injected current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(the_neuron, 1:11) = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate(k,1) = n(m(k));
end
R = rate(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('right_neurons.avi');
open(vidfile);
for i = 1:(length(Time)-1)
    % a movie of the model:
    plot(R/max(R));
    xlabel('neurons');
    ylabel('Rate');
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W_r*R));
    rate(:,i+1) = R;
end
close(vidfile);
hold off;


% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;





%%   QUESTION 3   %% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial injected current in [nA]:
I_inj_c = zeros(neurons_num,length(Time));
I_inj_l = zeros(neurons_num,length(Time));
I_inj_l(the_neuron, 200:300) = 100;
I_inj__r = zeros(neurons_num,length(Time));
I_inj_r(the_neuron, 600:800) = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate_c = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_c(k,1) = n(m(k));
end
R_c = rate_c(:,1);


% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate_l = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_l(k,1) = n(m(k));
end
R_l = rate_l(:,1);


% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate_r = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_r(k,1) = n(m(k));
end
R_r = rate_r(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('all_neurons.avi');
open(vidfile);
for i = 1:(length(Time)-1)
    R_c = R_c + dt*(1/Tau)*(-R_c + subplus(I_inj_c(:,i) + W*R_c + W*R_l + W*R_r));
    rate_c(:,i+1) = R_c;
    R_l = R_l + dt*(1/Tau)*(-R_l + subplus(I_inj_l(:,i) + W_l*R_l + W*rate_c(:,i) + W*R_r));
    rate_l(:,i+1) = R_l;
    R_r = R_r + dt*(1/Tau)*(-R_r + subplus(I_inj_r(:,i) + W_r*R_r + W*rate_c(:,i) + W*rate_l(:,i)));
    rate_r(:,i+1) = R_r;
    % a movie of the model:
    subplot(3,1,1);
    plot(R_c/max(R_c));
    xlabel('neurons');
    ylabel('Rate');
    subplot(3,1,2);
    plot(R_l/max(R_l));
    xlabel('neurons');
    ylabel('Rate');
    subplot(3,1,3);
    plot(R_r/max(R_r));
    xlabel('neurons');
    ylabel('Rate');
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
end
close(vidfile);
hold off;


%% plotting the activity of the chosen neuron VS time:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
plot(Time,rate(the_neuron,:));
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure', 'starting with different random firing rates'});
% ylim([-1 50]);
% lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
% title(lgd, 'neuron number:');
hold off;










%%%%%%%%%%%%%%%%%% what happens if one of the neurons gets an injected current and the initial rate of all neurons is zero? %%%%%%%%%%%%%%%%



% initial injected current in [nA]:
I_inj = zeros(neurons_num,length(Time));
I_inj(1, 1:10) = 40;


% initial rate of zero for all neurons in [Hz]:
rate = zeros(neurons_num,length(Time));
R = rate(:,1);


% the "rate" step (using the subplus function):
for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
end


% plotting the results until 50 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{13} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure (all starting at rest)', ['when only one of them (neuron 1) gets an injected current of ' num2str(I_inj(1,1)) 'nA for ' num2str(10*dt) 'sec (from the beginning)']});
ylim([-1 50]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;


% plotting the results until 100 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{13} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure (all starting at rest)', ['when only one of them (neuron 1) gets an injected current of ' num2str(I_inj(1,1)) 'nA for ' num2str(10*dt) 'sec (from the beginning)']});
ylim([-1 100]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;


% plotting the results until 150 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{13} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other in a ring structure (all starting at rest)', ['when only one of them (neuron 1) gets an injected current of ' num2str(I_inj(1,1)) 'nA for ' num2str(10*dt) 'sec (from the beginning)']});
ylim([-1 150]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;





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
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other (without any inhibition) in a ring structure', 'starting with different random firing rates'});
ylim([-1 50]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;


% plotting the results until 100 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other (without any inhibition) in a ring structure', 'starting with different random firing rates'});
ylim([-1 100]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;


% plotting the results until 150 Hz:
figure('Units','normalized','position',[0 0 1 1]);
hold on;
for i = 1:neurons_num 
    plot(Time,rate(i,:));
end
xlabel('Time [sec]');
ylabel('Rate [Hz]');
title({'\fontsize{14} Rate vs. Time', '\fontsize{12} \rm of 50 neurons affecting each other (without any inhibition) in a ring structure', 'starting with different random firing rates'});
ylim([-1 150]);
lgd = legend(string(1:neurons_num), 'Location', 'bestoutside', NumColumns = 2);
title(lgd, 'neuron number:');
hold off;





