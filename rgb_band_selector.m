function rgb_band_selector()
    %proper thread creation for the parfor loop
    clearvars
    disp(gpuDevice())

    if isempty(gcp('nocreate'))
        parpool('threads');
    end

    %this grabs the data files from the datacube
    %data_cube_og = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\data_cube.dat', [10880, 8160, 11], 'single', 0, 'bsq', 'ieee-le');        
    data_cube_og = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\vases_data_cube', [3625, 5449, 16], 'single', 0, 'bsq', 'ieee-le');        
    %data_cube_og = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\ege30_data_cube', [3397, 2389, 16], 'single', 0, 'bsq', 'ieee-le');        
    
    true_color = data_cube_og(:, :, [11, 7, 4]);
    lab_tci = rgb2lab(true_color);
    true_color_luma = lab_tci(:, :, 1);

    [pca_cube, ~, ~, ~] = pca_ims(data_cube_og);
    [res_log, num_components] = nosp();
    data_cube_og = pca_cube;%(:,:,1:num_components);

    for i=1:6
        figure;
        imagesc(pca_cube(:,:,i));
        colorbar;
        title(['Band ' num2str(i)]);
    end

    [~,~, bands] = size(data_cube_og);

    %downsampling for faster execution 8 is not a set number
    downsample_factor = 8;
    data_cube = imresize3(data_cube_og, [size(data_cube_og,1)/downsample_factor, size(data_cube_og,2)/downsample_factor, size(data_cube_og,3)], 'box');
    data_cube = im2double(data_cube);
    
    %calculates the array execution on one line like a lambda statement
    H = arrayfun(@(k) entropy(data_cube(:,:,k)), 1:bands);
    %disp('yurt');

    %this value can change
    bin_count = 100;

    %empty init
    mutual_information = zeros(bands, bands);

    %pre computing the edges and non bivariate histograms to save computation
    edges_all = cell(1,bands);
    hist_all = cell(1,bands);
    for i = 1:bands
        Bi = data_cube(:,:,i);
        edges_all{i} = linspace(min(Bi(:)), max(Bi(:)), bin_count+1);
        hist_all{i} = histcounts(Bi, edges_all{i}, 'Normalization', 'probability');
    end
    %disp('yurt');

    %threaded for loop, similar to parallel for each in c# but less work
    parfor i = 1:bands
        Bi = data_cube(:,:,i);
        ei = edges_all{i};
        histi = hist_all{i};

        %this is useful for threading, allows ease of access and elimiates some threading issues due to the indexing part
        local_mi = zeros(1,bands);
        
        for j = i+1:bands %no double calculations, mi is symmetric
            
            Bj = data_cube(:,:,j);
            ej = edges_all{j};
            histj = hist_all{j};
            
            %bivariate histogram. this sucker cant be lobbed onto the gpu and causes the need for parloops and downsampling
            histij = histcounts2(Bi, Bj, ei, ej, 'Normalization', 'probability');
            
            %some speedup found here from vectorization. the mask is so no 0% probabilites are run as a filtering mechanism
            mask = histij > 0;
            %these are necessary because the geniuses they employ over at mathworks decided it was a good idea to
            %collapse 'constant' matrices into one channel, so a 100x100 it deems similar can become a 100x1 or 1x100
            %there might be a simpler way to do this, but I am no expert
            pi_vec = histi(:); 
            pj_vec = histj(:);
            [i_idx, j_idx] = find(mask);

            %vectorization of a for loop with the use of mask and index checks to avoid matrix mulitplication issues
            MI = sum(histij(mask) .* log(histij(mask) ./ (pi_vec(i_idx) .* pj_vec(j_idx)))); %sig p(i,j) * log(p(i,j)/p(i)*p(j))
            
            %saving here and 'copying' mirror value later
            local_mi(j) = MI;
        end
        %good for threading to assign full values afterwards like this
        mutual_information(i,:) = local_mi;
    end

    %cool way to add the other half of data. since you full in one 'side' ~ its a triangle ~ the transformation properly fills in the other side neat.
    mutual_information = mutual_information + mutual_information';
    %disp('yurt');
    scores = H - sum(mutual_information,2)'/(bands-1); %saving a for loops f = H - 1/N-1 * sig (MI)
    scores = gather(scores);  
    disp(scores);

    
    %band_variance = band_variance(1:num_components);
    %scores = scores(1:num_components);

    %score_under = scores .* band_variance;
    %score_under = score_under/max(score_under);
    %score_over = scores ./ band_variance;
    %score_over = score_over/max(score_over);

    %band = sum(band_variance)/length(band_variance);
    %mean_closest = abs(band_variance - band);
    %true_color = mat2gray(rgb2gray(cat(3, data_cube_og(:,:,7), data_cube_og(:,:,5), data_cube_og(:,:,2))));
    %figure;
    %imshow(true_color);
    %var_true = sum(var(true_color,0, 'all'));
    %idea = abs(band_variance-var_true);
    %disp(idea);
    %mean_closest = (mean_closest/max(mean_closest));

    %this here is just testing to try and find the columns used to create the false color image for what I am working with
    %[mean_values, mean_order] = sort(mean_closest, 'descend');
    %[under_values, under_order] = sort(score_under,'descend');
    %[over_values, over_order] = sort(score_over,'descend');

    %disp(under_values);
    %disp(mean_values);
    %disp(over_values);
    %disp(under_order);
    %disp(over_order);
    %disp(mean_order);

    %{

    potential_bands = zeros(6,1);
    for i=1:3
        potential_bands((i-1)*2+1) = under_order(i);
        %potential_bands((i-1)*2+1) = over_order(i);
        potential_bands((i-1)*2+2) = mean_order(i);
    end

    potential_bands = unique(potential_bands);

    
    if under_values(1)/under_values(2) >= 2
        potential_bands = [potential_bands; under_order(1)];
    end

    %{
    if over_values(1)/over_values(2) >= 2
        potential_bands = [potential_bands; over_order(1)];
    end
    %}

    if mean_values(1)/mean_values(2) >= 2
        potential_bands = [potential_bands; mean_order(1)];
    end
    %disp(potential_bands);
    %}
    potential_bands = 1:num_components;
    combos = nchoosek(potential_bands,3);
    %nCombos = size(combos, 1);

    %disp('yurt');
    %potential_bands_gpu = gpuArray(potential_bands);
    data_cube_gpu = gpuArray(data_cube_og);
    %true_color = gpuArray(cat(3, data_cube_gpu(:,:,7), data_cube_gpu(:,:,5), data_cube_gpu(:,:,2)));
    
    values = zeros(1, 6);
    for i=1:6
        image = data_cube_gpu(:,:,i);
        image = abs(image);
        image = (image-min(image(:)))/(max(image(:))-min(image(:)));

        figure;
        imagesc(image);
        colorbar;

        %T = prctile(abs(image), 95);
        %mask = abs(image) > T;

        %figure;
        %imagesc(mask);
        %colorbar;
        %diff = var(double(res_log(mask)-true_color_luma(mask)), 0, 'all');
        res_mask = res_log .* image;
        tci_mask = true_color_luma .* image;
        output_image = res_mask - tci_mask;
        diff = var(double(output_image), 0, 'all');

        figure;
        imagesc(output_image);
        colorbar;
        
        values(i) = diff;
    end
    disp('here')
    disp(values)
    v = values .* scores;
    disp(v);

    differences = arrayfun(@(r,g,b) compute_diff(r,g,b, data_cube_gpu, res_log), combos(:,1), combos(:,2), combos(:,3));
    
    %disp('yurt');
    %ascend for underdrawing, descend for else
    [values, order] = sort(differences, 'ascend');
    output = combos(order, :);
    disp(values);
    disp(output);
    %min for underdrawing, max for else
    [~, idx_max] = min(differences);
    %[~, idx_max] = min(differences);
    %mean_val = mean(differences);
    %[~, idx_max] = min(abs(differences - mean_val));
    best_combo = combos(idx_max, :);
    best_r = best_combo(1);
    best_g = best_combo(2);
    best_b = best_combo(3);

    %disp('yurt');
    best_images = {
        cat(3, data_cube_gpu(:,:,best_r), data_cube_gpu(:,:,best_g), data_cube_gpu(:,:,best_b));
        cat(3, data_cube_gpu(:,:,best_r), data_cube_gpu(:,:,best_b), data_cube_gpu(:,:,best_g));
        cat(3, data_cube_gpu(:,:,best_g), data_cube_gpu(:,:,best_b), data_cube_gpu(:,:,best_r));
        cat(3, data_cube_gpu(:,:,best_g), data_cube_gpu(:,:,best_r), data_cube_gpu(:,:,best_b));
        cat(3, data_cube_gpu(:,:,best_b), data_cube_gpu(:,:,best_g), data_cube_gpu(:,:,best_r));
        cat(3, data_cube_gpu(:,:,best_b), data_cube_gpu(:,:,best_r), data_cube_gpu(:,:,best_g))
    };

    figure;
    for i = 1:6
        image = gather(best_images{i});
        subplot(2, 3, i);  
        imshow(image);
        title(['Image ' num2str(i)]);
    end

    figure;
    image = best_images{1};
    imshow(image);
    %save_curr_fig_to_file('BEST_BANDS_all_possible_outputs.jpg');
    disp(best_r);
    disp(best_g);
    disp(best_b);

    image = image.^.2;
    figure;
    imshow(image)


    
    %additional work required
    % it need some type of color science contrast color abs dif calculation when it picks the possible connections, and need 
    % a way to find out when to use multiple of the same channel. More information is required, especially with the color science 
    % stuff. But I think MI is a good road forward for band filtering. orthoganl subspace projection

end