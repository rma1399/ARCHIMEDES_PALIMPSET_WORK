function k_means_archimedes_palimpset_digital(img, k)
    
    %setting the image up for the best clustering use lab and high contrast
    %highlights
    im_double = im2double(img);
    im_contrast = im_double .^ 0.6;               % enhance contrast
    im_lab = rgb2lab(im_contrast);

    [height, width, ~] = size(im_lab);
    pixels = reshape(im_lab, height*width, 3);
    pixels = single(pixels);                      

    opts = statset('MaxIter', 10000, 'UseParallel', true);
    idx = kmeans(pixels, k, 'Options', opts, 'Replicates', 3, 'Start', 'plus');

    for i = 1:k
        filt = reshape(idx == i, height, width);   % logical mask
        filt_3d = cat(3, filt, filt, filt);       % replicate mask for RGB

        cluster_img = im_lab .* filt_3d;          % keep only cluster pixels
        max_val = max(cluster_img(:));
        if max_val > 0
            cluster_img = cluster_img / max_val;  % normalize for display
        end

        figure;
        imshow(cluster_img);
        title(['Cluster ', num2str(i)]);
    end
end