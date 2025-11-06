function [pca_output, num_components, eig_vectors, mu] = pca_ims(data_cube)
    %disp(gpuDevice());

    data_cube = gpuArray(data_cube);
    [rows, cols, bands] = size(data_cube);

    data = reshape(data_cube, [], bands);

    mu = mean(data,1);
    centered_data = data-mu;

    centered_data_covariance = (centered_data' * centered_data) / (size(centered_data, 1) - 1);

    [e, d] = eig(centered_data_covariance, 'vector');
    [eig_vals, order] = sort(d, 'descend');

    e = e(:, order);
    %disp(eig_vals);
    explained = eig_vals/sum(eig_vals) * 100;
    
    cum_var = cumsum(gather(explained)); %gather is used to %display my value for testing
    %disp(cum_var);

    num_components = find(cum_var >= 98, 1, 'first');
    %disp(num_components);

    batch_size = 1e5;
    num_batches = ceil(size(centered_data,1)/batch_size);
    score = zeros(size(centered_data,1), num_components, 'single', 'gpuArray');

    e = e(:, 1:num_components);

    for i = 1:num_batches
        idx = (i-1)*batch_size + 1 : min(i*batch_size, size(centered_data,1));
        score(idx,:) = centered_data(idx,:) * e;
    end

    pca_output = reshape(score, rows, cols, num_components);

    pca_output = gather(pca_output);
    eig_vectors = gather(e);
    mu = gather(mu);


    clear data_cube data centered_data centered_data_covariance d eig_vals order explained cum_var score idx

end