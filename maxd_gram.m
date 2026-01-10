function [endmembers, endmembers_index, volume] = maxd_gram(data, num, mnf_data, gram)
% data        : [npixels x nbands]
% num         : number of endmembers
% mnf_data    : MNF-transformed data (or 0 to reuse data)
% gram        : 'local' or 'general'

fprintf('---> In MaxD extracting endmembers and Grammian ...\n');

if isequal(mnf_data, 0)
    mnf_data = data;
end

% Transpose to match Python behavior
data  = data.';       % [nbands x npixels]
data2 = mnf_data.';   % [nbands x npixels]

[num_bands, num_pix] = size(data);

%% Magnitude of vectors
magnitude = sum(data.^2, 1);

idx1 = find(magnitude == max(magnitude));
idx2 = find(magnitude == min(magnitude));

% If multiple maxima/minima, keep first (matches Python behavior later)
idx1 = idx1(1);
idx2 = idx2(1);

%% Output arrays
endmembers       = zeros(num_bands, num);
endmembers_index = zeros(1, num);

% First two endmembers
endmembers(:,1) = data(:, idx1);
endmembers(:,2) = data(:, idx2);

endmembers_index(1) = idx1;
endmembers_index(2) = idx2;

%% Projection setup
data_proj = data2;
I = eye(num_bands);

volume = zeros(1, num);

%% Main loop
for i = 3:num

    % Difference vector
    diff = data_proj(:, idx2) - data_proj(:, idx1);   % [nbands x 1]

    % Pseudoinverse
    pseudo = pinv(diff);

    % Orthogonal projection
    data_proj = (I - diff * pseudo) * data_proj;

    % Update reference index
    idx1 = idx2;

    % Distance computation
    diff_new = sum(((data_proj(:, idx2) * ones(1, num_pix)) - data_proj).^2, 1);

    % Find next endmember
    idx2 = find(diff_new == max(diff_new));

    % Handle multiple maxima
    if numel(idx2) > 1
        idx2 = idx2(1);
    end

    % Assign endmember
    endmembers(:, i-1) = data(:, idx2);
    endmembers_index(i-1) = idx2;

    if strcmp(gram, 'general')
        gen_gram = endmembers(:,1:i-1).' * endmembers(:,1:i-1);
        volume(i-1) = det(gen_gram); %sqrt(abs(det(gen_gram)));
    end
end
end
