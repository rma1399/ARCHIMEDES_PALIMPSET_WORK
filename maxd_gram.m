function [endmembers, endmembers_index, volume] = maxd_gram(data, num, mnf_data, gram)
% data        : [npixels x nbands]
% num         : number of endmembers
% mnf_data    : MNF-transformed data (or 0 to reuse data)
% gram        : 'local' or 'general'

fprintf('---> In MaxD extracting endmembers and Grammian ...\n');

if isequal(mnf_data, 0)
    mnf_data = data;
end

data  = data.';       
data2 = mnf_data.';   

[num_bands, num_pix] = size(data);

magnitude = sum(data.^2, 1);

idx1 = find(magnitude == max(magnitude));
idx2 = find(magnitude == min(magnitude));

idx1 = idx1(1);
idx2 = idx2(1);

endmembers       = zeros(num_bands, num);
endmembers_index = zeros(1, num);

endmembers(:,1) = data(:, idx1);
endmembers(:,2) = data(:, idx2);

endmembers_index(1) = idx1;
endmembers_index(2) = idx2;

data_proj = data2;
I = eye(num_bands);

volume = zeros(1, num);

for i = 3:num

    diff = data_proj(:, idx2) - data_proj(:, idx1); 
    pseudo = pinv(diff);

    data_proj = (I - diff * pseudo) * data_proj;

    idx1 = idx2;
    diff_new = sum(((data_proj(:, idx2) * ones(1, num_pix)) - data_proj).^2, 1);
    idx2 = find(diff_new == max(diff_new));
    if numel(idx2) > 1
        idx2 = idx2(1);
    end
    
    endmembers(:, i-1) = data(:, idx2);
    endmembers_index(i-1) = idx2;

    if strcmp(gram, 'general')
        gen_gram = endmembers(:,1:i-1).' * endmembers(:,1:i-1);
        volume(i-1) = det(gen_gram); %sqrt(abs(det(gen_gram)));
    end
end
end
