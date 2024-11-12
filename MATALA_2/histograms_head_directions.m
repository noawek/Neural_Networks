


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


         head_direction_degrees{k} = head_direction_deg;
         head_direction_degrees{k+1} = head_direction_deg;


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


    head_direction_degrees{k} = head_direction_deg;


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




figure('Units','normalized','position',[0 0 1 1]);
hold on;    
subplot(3,2,1);
histogram(head_direction_degrees{1});
title('neuron 1');
subplot(3,2,2);
histogram(head_direction_degrees{2});
title('neuron 2');
subplot(3,2,3);
histogram(head_direction_degrees{3});
title('neuron 3');
subplot(3,2,4);
histogram(head_direction_degrees{4});
title('neuron 4');
subplot(3,2,5);
histogram(head_direction_degrees{5});
title('neuron 5');
subplot(3,2,6);
histogram(head_direction_degrees{6});
title('neuron 6');
sgtitle('Head Directions while measuring each neuron:')
hold off;