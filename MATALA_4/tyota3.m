



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



% things related to time in [sec]:
Tau = 10;  
dt = 0.1;
Time = 0:dt:100;





%%   QUESTION 2 (A) - creating a ring with a diversion to the left  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_l = zeros(neurons_num);


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_l,1)
    for j = 1:size(W_l,2)
        Di_j = min(abs((j-5)-i), neurons_num - abs((j-5)-i));
        W_l(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end




%%   QUESTION 2 (B) - creating a ring with a diversion to the right  %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_r = zeros(neurons_num);


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_r,1)
    for j = 1:size(W_r,2)
        Di_j = min(abs((j+5)-i), neurons_num - abs((j+5)-i));
        W_r(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
    end
end



%%   QUESTION 3   %% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial injected current in [nA]:
I_inj_l = zeros(neurons_num,length(Time));
I_inj_l(:, 100:200) = 10;
I_inj_r = zeros(neurons_num,length(Time));
I_inj_r(:, 400:500) = 10;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the initial random conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initial random rates of neurons in [Hz]:
m = randperm(neurons_num);
n = 0.025:0.025:5;
rate_c = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_c(k,1) = n(m(k));
end
rate_c(:,1) = rate_c(:,1)/max(rate_c(:,1));
R_c = rate_c(:,1);
rate_l = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_l(k,1) = n(m(k));
end
rate_l(:,1) = rate_l(:,1)/max(rate_l(:,1));
R_l = rate_l(:,1);
rate_r = zeros(neurons_num,length(Time));
for k = 1:neurons_num
    rate_r(k,1) = n(m(k));
end
rate_r(:,1) = rate_r(:,1)/max(rate_r(:,1));
R_r = rate_r(:,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the calculation and representation of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the "rate" step (using the subplus function):
figure('Units','normalized','position',[0 0 1 1]);
vidfile = VideoWriter('allofthe_neurons.avi');
open(vidfile);
for i = 1:(length(Time)-1)
    % a movie of the model:
    subplot(3,1,1);
    plot(rate_c(:,i));
    xlabel('neurons');
    ylabel('Rate');
    subplot(3,1,2);
    plot(rate_l(:,i));
    ylim([min(rate_l(:,i)) 1]);
    xlabel('neurons');
    ylabel('Rate');
    subplot(3,1,3);
    plot(rate_r(:,i));
    ylim([min(rate_r(:,i)) 1]);
    xlabel('neurons');
    ylabel('Rate');
    sgtitle(['time = ' num2str(floor((i-1)*dt))]);
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
    % rate step:
    R_c = R_c + dt*(1/Tau)*(-R_c + subplus(W*R_c + W*R_l + W*R_r));
    R_c = R_c/max(R_c);
    rate_c(:,i+1) = R_c;
    R_l = R_l + dt*(1/Tau)*(-R_l + subplus(W_l*R_l + W*rate_c(:,i) + W*R_r) + I_inj_l(:,i));
    R_l = R_l/max(R_l);
    rate_l(:,i+1) = R_l;
    R_r = R_r + dt*(1/Tau)*(-R_r + subplus(W_r*R_r + W*rate_c(:,i) + W*rate_l(:,i)) + I_inj_r(:,i));
    R_r = R_r/max(R_r);
    rate_r(:,i+1) = R_r;
end
close(vidfile);
hold off;


