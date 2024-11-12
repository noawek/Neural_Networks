

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setteing the neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setting some variables:
W_c_l = zeros(neurons_num);
W_c_r = zeros(neurons_num);
W_r_l = zeros(neurons_num);


% the W's matrix (setting the sinapses strengths and directions between each pair of neurons and neurons to themselves):
for i = 1:size(W_r_l,1)
    for j = 1:size(W_r_l,2)
        Di_j = min(abs(j-i), neurons_num - abs(j-i));
        W_c_l(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
        W_c_r(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;
        W_r_l(i,j) = exp(-((Di_j^2)/(sigma1^2)))-0.1;        
    end
end

%%


% the "rate" step (using the subplus function):
for i = 1:(length(Time)-1)
    R = R + dt*(1/Tau)*(-R + subplus(I_inj(:,i) + W*R));
    rate(:,i+1) = R;
    % a movie of the model:
    plot(R/max(R));
    xlabel('neurons');
    ylabel('Rate');
    pause(0.05);
end