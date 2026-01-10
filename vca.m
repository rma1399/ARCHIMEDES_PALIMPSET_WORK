function S = vca(data, end_members)
    %finds pure endmembers to use for NOSP
    [~, bands] = size(data);

    data = data'; %transposes the data to make calculations easier
    %centered_data = data - mean(data, 2);
    centered_data = data;

    %this is to remove dimentionality and project the basis vectors
    [U, ~, ~] = svd(centered_data * centered_data', 'econ'); 
    Ud = U(:, 1:end_members-1);       
    Y = Ud' * centered_data;           
    
    disp('here');
    S = zeros(bands, end_members);
    for i = 1:end_members
        disp('here');
        f = randn(end_members-1, 1); %point of randomness inside vca which can make reproducibility difficult
        f = f / norm(f);  

        %Project pixels onto f
        projections = f' * Y;

        %Pick the extreme pixel
        [~, idx] = max(abs(projections)); %getting our pure endmembers
        S(:, i) = centered_data(:, idx); %assign idx values to our S output
    end
end