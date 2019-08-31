    clear all
    close all
    
    %% setup variables
    
    time_frequency = 32; %Hz
    Hz = 32;
    
    file_name = 'file_name';
    
    threshold = 15; %for event detection
    %In total_amplitude, it should be 0 otherwise, it will be 0 when it's below 10%
    
    % event display
    m = 30; % how seperate from each data
    m_2 = 1.01; % Where event dot is located (at the max amp or each ROI)
    markersize = 300;
    
    %% LOAD
      
    table(:,:) = readmatrix(file_name);
    f_signal=table'; %if trasposing is needed
    [r, c, z] = size(f_signal);
    sec_length = r/time_frequency;
    xval = 1/time_frequency:1/time_frequency:sec_length; % sec for x axis

%% figure; plot(f_signal)
    
    for i = 1:z
        figure
        plot(xval, f_signal(:,:,i))
        title(i); xlabel('Time (sec)');
        
        figure; imagesc(f_signal'.*100); colorbar; %intensity profile plotting
    end
    
 % f_signal display without overlapping
    f_display = zeros(r,c,z);
    for i = 1:z
        for j = 1:c
            f_display(:,j,i) = f_signal(:,j,i) -m*(j-1);
        end
    end

    for i = 1:z
        figure
        plot(xval, f_display(:,:,i))
        title(i); xlabel('Time (sec)');
    end

    %% correlationi matrix 
    % + average correlation 
    average_correlation = zeros(z,1);

    %check the minimum correlation values across the sheets (min_corr)
    min_corr = 1;
    max_corr = 0;
    for i = 1:z
        min_corr = min([min(min(corrcoef(f_signal(:,:,i)))), min_corr]);

        ascend_sort  = sort(corrcoef(f_signal(:,:,i)), 'descend');
        max_corr = max(max(ascend_sort(2,:)), max_corr);
     end

    figure
    for i = 1:z
        % cell by cell correlation
        cell_correlation = corrcoef(f_signal(:,:,i));

        subplot(2,z,i);
        imagesc(cell_correlation,[min_corr max_corr])
        title(i)
        average_correlation(i) = mean(mean(cell_correlation));


        % time by time correlation
        time_correlation = corrcoef(f_signal(:,:,i)');

        subplot(2,z,i+z);
        imagesc(time_correlation)
        title(i)
    end

    figure
    plot(average_correlation,'ro--');
    title('average correlation')

    %% event detection
    %%%%%%%%%%%%%%% make raster plot + PSTH
    raster_matrix = zeros(size(f_signal));
    raster_matrix_pre = zeros(size(f_signal));
    [r c z] = size(f_signal);
   
    
    for i = 1:z
        for j = 1:c
            [q w] = findpeaks(f_signal(:,j,i));
            raster_matrix_pre(w,j,i) = q;
        end
    end
    
    temporal_code = zeros(z,r);
    raster_matrix(raster_matrix_pre > threshold) = 1;
    
    raster_matrix_amplitude = raster_matrix_pre;
    raster_matrix_amplitude(raster_matrix_pre <= threshold) = 0;
    
    for i = 1:z;
        figure;
        subplot(2,1,1);
        imshow(not(raster_matrix(:,:,i)'));
        title(i);

        subplot(2,1,2);
        temporal_code(i,:) = sum(raster_matrix(:,:,i)');
        bar(temporal_code(i,:)); 
        axis tight;
        axis([min(xval) max(xval*32) 0 max(temporal_code)]);
    end;

   
    %% event disply
    event_display_matrix = zeros(r,c*2,z);
    
    
    for j = 1:z
        for i = 1:2:c*2
            event_display_matrix(:,i+1,j) = f_signal(:,round(i/2),j);
            flag_position = max(event_display_matrix(:,i+1,j))*m_2;
            event_display_matrix(:,i,j) = raster_matrix(:,round(i/2),j)*flag_position;
        end
    end
    
    %% Use display function
    event_display(event_display_matrix, [1:c], m, xval, markersize);
    % ROI index can be assinged to the second variable %[1:C]==> all roi
    
    %%  ROI ? rate frequency
    
    rate_code = zeros(c,N_compartment, z); %c = ROI number
            
    for i = 1:z
        for j = 1:N_compartment
            x = find(xval == compart_info(j,1));
            y = find(xval == compart_info(j,2));
            rate_code(:,j,i) = sum(raster_matrix(x:y,:,i))/(compart_info(j,2)-compart_info(j,1));
        end
    end
    
    
    xy_scale_ratecode = [0 c+1 0 max(max(rate_code))+0.5]; %column?? [Xmin Xmax Ymin Ymax]
    
    for i = 1:z
        for j = 1:N_compartment

        figure
        bar(rate_code(:,j,i))
        axis(xy_scale_ratecode);
        title(['Event frequency ' num2str(i) '  -  ' num2str(j)]) 
        end
        
    end
    
    %% ROI ? total amplitude (??? %, event ??? ?? ?? amplitude?)
    total_amplitude = zeros(c,N_compartment, z);
    
for i = 1:z
        for j = 1:N_compartment
            x = find(xval == compart_info(j,1));
            y = find(xval == compart_info(j,2));
            total_amplitude(:,j,i) = sum(raster_matrix_amplitude(x:y,:,i));
            num = sum(raster_matrix_amplitude(x:y,:,i)>0)';
            total_amplitude(:,j,i) = total_amplitude(:,j,i)./num;
            total_amplitude(:,j,i) = total_amplitude(:,j,i) .* 100;
            
        end
    end

    
    xy_scale_totalamplitude = [0 c+1 0 max(max(total_amplitude))+3]; %column?? [Xmin Xmax Ymin Ymax]
    
    for i = 1:z
        for j = 1:N_compartment

        figure
        bar(total_amplitude(:,j,i))
        axis(xy_scale_totalamplitude);
        title(['Total amplitudes ' num2str(i) '  -  ' num2str(j)]) 
        end
        
    end    
    
    
%% histogram? ?? ?? (normalized rate code / total amplitude? subtraction)

rate_mean = mean(rate_code,2);
rate_mean = repmat(rate_mean, 1,c);
rate_code_nor = rate_code./rate_mean;

amplitude_mean = mean(total_amplitude,2);
amplitude_mean = repmat(amplitude_mean, 1,c);
amplitude_code_nor = total_amplitude./amplitude_mean;

for i = 1:z-1
    figure
    bar(rate_code_nor(i+1,:) - rate_code_nor(1,:));
    title(['Ralative changes to spontaneous activity (event frequency) ' num2str(i+1)])
end
    
for i = 1:z-1
    figure
    bar(amplitude_code_nor(i+1,:) - amplitude_code_nor(1,:));
    title(['Ralative changes to spontaneous activity (total amplitude) ' num2str(i+1)])
end
    
%% ROI ?  display
roi_no = input('ROI list?\n');
cell_display(roi_no, r, z, m, f_signal);



%% SRoh
% csvwrite('R4_2_Image10_spon5min_CF_amp.csv',raster_matrix_amplitude)







