function diff = compute_diff(r, g, b, data_cube, res_log)

    % select the three bands
    img_r = data_cube(:,:,r);
    img_g = data_cube(:,:,g);
    img_b = data_cube(:,:,b);

    % create composite image
    image = cat(3, img_r, img_g, img_b);


    % difference to true color
    %diff = var(double(image), 0, 'all');
    %image_lab = rgb2lab(image);
    %flat = reshape(image, [], 3);
    %mu = mean(flat);
    %diff = mean(sum(flat-mu).^2, 2);
    %diff = mean(log(1 + (flat(:) - mean(flat(:))).^2));
    %diff = var(image, 0, 'all')-var(double(image-res_log),0,'all');
    res_im = sqrt(sum(image.^2, 3)); 
    diff = var(double(res_log - res_im), 0, 'all');
end