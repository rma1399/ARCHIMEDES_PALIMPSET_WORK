function diff = compute_diff(r, g, b, data_cube)

    % select the three bands
    img_r = data_cube(:,:,r);
    img_g = data_cube(:,:,g);
    img_b = data_cube(:,:,b);

    % create composite image
    image = cat(3, img_r, img_g, img_b);


    % difference to true color
    %diff = var(double(image), 0, 'all');
    %image_lab = rgb2lab(image);
    flat = reshape(image, [], 3);
    %mu = mean(flat);
    %diff = mean(sum(flat-mu).^2, 2);
    diff = mean(log(1 + (flat(:) - mean(flat(:))).^2));
end