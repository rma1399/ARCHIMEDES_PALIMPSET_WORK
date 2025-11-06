function MNF_output = mnf_ims(data_cube)
    %document this out with the math formulas more
    % Z = S + N
    disp(gpuDevice());
    data_cube = gpuArray(data_cube);
    
    [rows, cols, bands] = size(data_cube);

    data = reshape(data_cube, rows*cols, bands);

    mu = mean(data,1);
    centered_data = data-mu;

    %horizontal pixel difference, shoulds matter under assumption noise is not correlated to direction
    noise = centered_data(2:end, :) - centered_data(1:end-1,:); 
    noise = noise-mean(noise,1);


    batch_size = 1e5; 
    [num_rows, num_bands] = size(noise);
    cov_matrix = gpuArray(zeros(num_bands, num_bands, 'single')); 

    for i = 1:ceil(num_rows/batch_size)
        idx = (i-1)*batch_size + 1 : min(i*batch_size, num_rows);
        batch = noise(idx, :);
        cov_matrix = cov_matrix + batch' * batch;  % accumulate
    end
    noise_cov_matrix = cov_matrix / (num_rows - 1);



    [e, d] = eig(noise_cov_matrix);
    whitening_vector = e*diag(1 ./ sqrt(diag(d))) * e'; %<-- W = E*D^(-1/2)*E^T
    whitened_data = centered_data * whitening_vector'; % A*W^T
    whitened_data = whitened_data - mean(whitened_data, 1);


    w_cov_matrix = gpuArray(zeros(num_bands, num_bands, 'single'));
    for i = 1:ceil(num_rows/batch_size)
        idx = (i-1)*batch_size + 1 : min(i*batch_size, num_rows);
        batch = whitened_data(idx, :);
        w_cov_matrix = w_cov_matrix + batch' * batch;  % accumulate
    end
    
    whitened_cov_matrix = w_cov_matrix / (num_rows - 1);


    [e_MNF, d_MNF] = eig(whitened_cov_matrix);

    [eig_vals, order] = sort(diag(d_MNF), 'descend');
    e_MNF = e_MNF(:, order);

    explained = eig_vals/sum(eig_vals) * 100;
    cum_var = cumsum(gather(explained));

    disp(cum_var);
    num_components = bands; %find(cum_var >= 98, 1, 'first');

    disp(num_components);

    MNF = whitened_data * e_MNF(:,1:num_components); %order important because white_data(row*col, bands) and e_MNF is (bands, num_components)
    disp(size(MNF));

    MNF_output = reshape(MNF, rows, cols, num_components);

    MNF_output = gather(MNF_output);

    clear data_cube cov_matrix data centered_data noise noise_cov_matrix whitened_data whitened_cov_matrix whitening_vector e e_MNF explained eig_vals d d_MNF order cum_var
end