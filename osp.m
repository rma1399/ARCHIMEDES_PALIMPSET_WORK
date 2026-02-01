function band_variance = osp()

    %create the subspace with the datacube
    %evaluate the projection
    %label data such that you can identify the background
    %idenitfy unique orthoganal points
    %get rid of noisy unqiues and label as undertext
    %bonus: find which unique projection contains the undertext and label the data points
    
    %this grabs the data files from the datacube
    %input= imhypercube('data_cubes/ege30_data_cube.hdr');

    %gathers the data into it in real size
    %data_cube = gather(input);
    
    %info = enviinfo('data_cubes\ege30_data_cube.hdr');

    % Read the ENVI cube
    %data_cube = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\data_cube.dat', [10880, 8160, 11], 'single', 0, 'bsq', 'ieee-le');        
    data_cube = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\vases_data_cube', [3625, 5449, 16], 'single', 0, 'bsq', 'ieee-le');        
    %data_cube_og = multibandread('C:\Users\ryana\OneDrive\Documents\MATLAB\ARCHIMEDES_PALIMPSET_WORK\data_cubes\ege30_data_cube', [3397, 2389, 16], 'single', 0, 'bsq', 'ieee-le');        


    %{
    x1 = 200.00;
    y1 = 200.00;
    x2 = 2200.00;
    y2 = 3200.00;
    

    data_cube = data_cube_og(y1:y2, x1:x2, :);
    %}

    figure;
    imshow(data_cube(:,:,1));

    disp(size(data_cube));
    
    [pca_cube, num_components, eigen_vectors, mu] = pca_ims(data_cube);

    
    [rows, cols, bands] = size(data_cube);
    data = reshape(data_cube, [], bands);
    
    %data = data.^(0.4);
    centered_data = data-mu;
    Q = eigen_vectors(:,1:num_components);

    orthoganl_projection = eye(size(Q, 1)) - Q*Q'; % I - QQ^T

    projected_data = centered_data * orthoganl_projection';
    figure;
    plot(projected_data);

    osp_output = reshape(projected_data, rows, cols, bands);

    disp(num_components)
    figure;
    for i = 1:num_components
        subplot(1,num_components,i);
        imshow(mat2gray(osp_output(:,:,i)));
        title(['OSP Band ' num2str(i)]);
    end
    %save_curr_fig_to_file('OSP_OUTPUTS_1to3.jpg');
    

    osp_energy = sqrt(sum(osp_output.^2, 3));
    disp(size(osp_energy));
    figure;
    imagesc(osp_energy);
    axis image off;
    colormap viridis
    colorbar;
    title('OSP Residual Energy (Pigment Map)');
    save_curr_fig_to_file('OSP_vases_energy_map.jpg');

    band_variance = squeeze(var(reshape(osp_output, [], bands)));
    %band_variance = squeeze(var(reshape(pca_cube, [], bands)));
    %band_variance = squeeze(var(reshape(data_cube, [], bands)));
    figure;
    plot(band_variance);
    xlabel('Band'); ylabel('Residual variance');
    title('Bands sensitive to undertext');
    %save_curr_fig_to_file('UNDERTEXT_by_band.jpg');

    figure;
    image = cat(3, pca_cube(:,:,4), pca_cube(:,:,3), pca_cube(:,:,2));
    image = image.^0.2;
    imshow(image);

    %{
    [values, order] = sort(band_variance, "ascend");
    [~, bands] = size(projected_data);
    %find way to find different peaks or subtract the background
    for i=1:bands
        figure;
        h = histogram(projected_data(:,i)-projected_data(:, order(2)));
        disp(length(h.Values));
    end
    %}

    

    %{
    threshold = 0.05*osp_energy;
    pixel_data = osp_output(osp_energy>threshold);

    
    [idx,~] = kmeans(pixel_data, 3, 'MaxIter', 500, 'Replicates', 5);

    pigment_map = zeros(rows*cols, 1);
    pigment_map(osp_energy>threshold) = idx;

    pigment_map = reshape(pigment_map, rows, cols);
    figure;
    imagesc(pigment_map);
    axis image off;
    colormap(lines);  % distinct colors per pigment
    colorbar;
    title('Estimated pigment clusters');
    %}
    
    

end