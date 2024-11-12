



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% introduction to neuronal networks - assignment 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
    POSX{k} = posx;
    POSY{k} = posy;
    

    % setting the bins of locations for the analyses (each 5cmx5cm are grouped together):
    for i = 1:20
        small_squares(i,:) = [((i-1)*5), (i*5)];
    end
    
    
    % extracting some important variables from the files (using "if" to discriminate the first 4 files, each includes data only on one
    % neuron recorded, from the last file that includes data recorded from 2 neurons simultaniously):  

   
    % extracting spike trains of the two last neurons (in the last file), so we can calculate rates later on:
    if length(spiketrain_before) == length(spiketrain)   
    
        % loading some relevant data from them (about the location of the rat) before continuing:
        POSY{k+1} = posy;
        POSX{k+1} = posx;
        
        % neuron 5:
        spiketrains{k} = spiketrain1;
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
        spiketrains{k+1} = spiketrain2;
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
         

         % permuting the spike trains data as a preperation for Question 3: 

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
         for j = 1:length(small_squares)
             for i = 1:length(small_squares)
             places = find(POSY{k} < small_squares(j,2) & POSY{k} > small_squares(j,1) & POSX{k} < small_squares(i,2) & POSX{k} > small_squares(i,1));
             countperpos = length(places); 
             timeperpos = countperpos*(1/sampleRate);
             num_spikes1 = sum(spiketrain1(places));
             num_spikes2 = sum(spiketrain2(places));
             rateperpos(k,j,i) = num_spikes1/timeperpos;
             rateperpos(k+1,j,i) = num_spikes2/timeperpos;
             num_spikes1_fake = sum(spiketrain_fakes{k}(places));
             num_spikes2_fake = sum(spiketrain_fakes{k+1}(places));
             rateperpos_fake(k,j,i) = num_spikes1_fake/timeperpos;
             rateperpos_fake(k+1,j,i) = num_spikes2_fake/timeperpos;
             num_spikes1_fake2 = sum(spiketrain_fakes2{k}(places));
             num_spikes2_fake2 = sum(spiketrain_fakes2{k+1}(places));
             rateperpos_fake2(k,j,i) = num_spikes1_fake2/timeperpos;
             rateperpos_fake2(k+1,j,i) = num_spikes2_fake2/timeperpos;
             end
         end
    
         
    % extracting the same important variables and preforming relevant calculations, 
    % but this time for the rest of the neurons (1-4):
    else


    % zeroing the APs variable each the the loop runs (so data won't mix):
    APs = 0;
    
    % extracting spike trains of those neurons (from the first 4 files), so we can calculate rates later on:
    spiketrains{k} = spiketrain;
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


    % permuting the spike trains data as a preperation for Question 3:

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
    for j = 1:length(small_squares)
        for i = 1:length(small_squares)
        places = find(POSY{k} < small_squares(j,2) & POSY{k} > small_squares(j,1) & POSX{k} < small_squares(i,2) & POSX{k} > small_squares(i,1));
        countperpos = length(places); 
        timeperpos = countperpos*(1/sampleRate);
        num_spikes = sum(spiketrain(places));
        num_spikes_fake = sum(spiketrain_fakes{k}(places));
        num_spikes_fake2 = sum(spiketrain_fakes2{k}(places));
        rateperpos(k,j,i) = num_spikes/timeperpos;
        rateperpos_fake(k,j,i) = num_spikes_fake/timeperpos;
        rateperpos_fake2(k,j,i) = num_spikes_fake2/timeperpos;
        end
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



for i = 1:size(rateperpos,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    plot(POSX{i},POSY{i});  
    scatter(POSX{i}(find(spiketrains{i} == 1)),POSY{i}(find(spiketrains{i} == 1)),15,'red','filled');
    axis square;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    legend('the rat track', 'action potential','Location', 'bestoutside');
    title({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Superimposing the Spatial Activity of the Neuron [Action Potentials] on the Spatial Activity of the Rat [her track in a 100cmX100cm squared environment]'});
    hold off;
end


% clearing some irrelevant variables so the workspace would be comfetrable to work with:
clear num_spikes_fake2 num_spikes2_fake2 num_spikes1_fake2 AP_permute AP_permute1 AP_permute2 F AP_perm1 AP_perm2 APs APs1 APs2 APs2_before APs1_before APs_before boxSize countperpos headDirection sampleRate spiketrain1 spiketrain2 spiketrain_fake spiketrain AP_perm ans j b m k c l a n i num num_spikes2_fake num_spikes1_fake num_spikes num_spikes_fake num_spikes2 num_spikes1 spiketrain_before timeperpos posx posy post theta phase places;





%%   QUESTION 2   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HitMaps of all 6 neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:size(rateperpos,1)
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    imagesc([2.5 97.5], [2.5 97.5], imgaussfilt(squeeze(rateperpos(i,:,:)))); 
    xlim([0 100]);
    ylim([0 100]);
    axis square;
    colormap("hot");
    c = colorbar;
    % setting the same limits for the colors scale for all neurons, so we can compare them:
    clim([0 54]);         
    c.Label.String = 'Firing Rate [Hz]';
    c.Label.FontSize = 12;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rate [Hz] as a function of the Location of the rat [100cmX100cm squared environment]'});
    hold off;
end


% clearing some irrelevant variables so the workspace would be comfetrable to work with:
clear i c




%%   QUESTION 3   %%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures of all 6 neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comparing the "fake" (randomized) data to the original one %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% a small reminder: 
% we have already managed to get a vector of the margines between action
% potentials, played with it while permuting it twice, and created 2 new
% spike trains for each neuron according to the new order of margines dictated
% after that, we've calculated firing rates per position for each permutation... 


% after setting the bed for it while loading and re-ordering the data in the big loop of Question 1,
% now it is the time to examine and compare the randomized rates to the original rates of all 6 neurons:
for i = 1:size(rateperpos_fake,1)
    % a comparison between graphs similar to those shown in question 1 (just for fun.. it's not on the presentation file):
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Superimposing the Spatial Activity of the Neuron [Action Potentials] on the Spatial Activity of the Rat [her track in a 100cmX100cm squared environment]:', 'comparison between the original data from the neuron, and the test data (where the inter-spike interval of the spike train was randomly permuted)', '\fontsize{10} to check if a neuron can really be classified either as a place cell or as a grid cell'});
    subplot(1,3,1);
    hold on;
    plot(POSX{i},POSY{i});  
    scatter(POSX{i}(find(spiketrains{i} == 1)),POSY{i}(find(spiketrains{i} == 1)),10,'red','filled');
    axis square;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title('Original Data');
    subplot(1,3,2);
    hold on;
    plot(POSX{i},POSY{i});  
    scatter(POSX{i}(find(spiketrain_fakes{i} == 1)),POSY{i}(find(spiketrain_fakes{i} == 1)),10,'red','filled');
    axis square;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title('Randomly Permuted Data 1');
    subplot(1,3,3);
    hold on;
    plot(POSX{i},POSY{i});  
    scatter(POSX{i}(find(spiketrain_fakes2{i} == 1)),POSY{i}(find(spiketrain_fakes2{i} == 1)),10,'red','filled');
    axis square;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    % one legend (relevent for all 3 garphs):
    Lgdn = legend('the rat track', 'action potential','Location', 'bestoutside');
    Lgdn.Position(1) = 0.48;
    Lgdn.Position(2) = 0.8;
    title('Randomly Permuted Data 2');
    hold off;
    % a comparison between graphs similar to those shown in question 2:
    figure('Units','normalized','position',[0 0 1 1]);
    hold on;
    sgtitle({['\fontsize{14} \bf Neuron ' num2str(i) ' -'], '\fontsize{12} \rm Firing Rates [Hz] as a function of the Location of the rat [100cmX100cm squared environment]:', 'comparison between the original data from the neuron, and the test data (where the inter-spike interval of the spike train was randomly permuted)', '\fontsize{10} to check if a neuron can really be classified either as a place cell or as a grid cell'});
    subplot(1,3,1);
    hold on;
    imagesc([2.5 97.5], [2.5 97.5], imgaussfilt(squeeze(rateperpos(i,:,:)))); 
    xlim([0 100]);
    ylim([0 100]);
    yticks(0:20:100);
    axis square;
    % setting the same limits for the colors scale for all neurons, so we can compare them:
    colormap("hot");
    clim([0 54]);
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title('Original Data');
    subplot(1,3,2);
    hold on;
    imagesc([2.5 97.5], [2.5 97.5], imgaussfilt(squeeze(rateperpos_fake(i,:,:)))); 
    xlim([0 100]);
    ylim([0 100]);
    yticks(0:20:100);
    axis square;
    % setting the same limits for the colors scale for all neurons, so we can compare them:
    colormap("hot");
    clim([0 54]);
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title('Randomly Permuted Data 1');
    subplot(1,3,3);
    hold on;
    imagesc([2.5 97.5], [2.5 97.5], imgaussfilt(squeeze(rateperpos_fake2(i,:,:)))); 
    xlim([0 100]);
    ylim([0 100]);
    yticks(0:20:100);
    axis square;
    colormap("hot");
    % one colorbar (relevent for all 3 garphs):
    c = colorbar;
    clim([0 54]);     % setting the same limits for the colors scale for all neurons, so we can compare them:
    c.Position(1) = 0.94;
    c.Position(2) = 0.35;
    c.Label.String = 'Firing Rate [Hz]';
    c.Label.FontSize = 12;
    ylabel('position on the y-axis [cm]');
    xlabel('position on the x-axis [cm]');
    title('Randomly Permuted Data 2');
    hold off;
end


% clearing some irrelevant variables so the workspace would be comfetrable to work with:
clear i








