function make_data_cube()
    reset(gpuDevice());
    disp(gpuDevice());
    %%images
    
    disp('not_cache')
    processed_images = {
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED365_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED445_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED470_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED505_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED530_01_raw.tif')));
        mat2gray(rgb2gray(imread("081r-088v_Arch03r_Sinar_LED570_01_raw.tif")));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED617_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED625_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED700_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED735_01_raw.tif')));
        mat2gray(rgb2gray(imread('081r-088v_Arch03r_Sinar_LED870_01_raw.tif')))
    };
    disp('test')


    %%data cube
    disp('making data cube ...');
    num_of_bands = length(processed_images);
    for i=1:num_of_bands
        processed_images{i} = im2single(processed_images{i});
    end

    pcd = cat(3, processed_images{:});
    data_cube = gpuArray(pcd);
    %true_color = gpuArray(cat(3, processed_images{7}, processed_images{5}, processed_images{2}));
    clear processed_images
    disp('done');
    

    %%pca
    disp('making pca images ...');
    %pca_cube = pca_ims(data_cube);
    disp('done');
    
    

    %%mnf
    disp('making mnf images ...');
    %mnf_cube = mnf_ims(data_cube);
    disp('done');
    
    %%make cubes into list
    [rows, cols, bands] = size(data_cube);
    full_cube = zeros(rows, cols, bands,'single');
    
    full_cube(:,:,1:bands) = data_cube;
    %full_cube(:,:,1+bands:2*bands) = pca_cube;
    %full_cube(:,:,1+bands*2:end) = mnf_cube;


    enviwrite(imhypercube(full_cube), 'data_cube.envi');
    reset(gpuDevice());

    %%my work current
    %rgb_bands_differencer(full_cube, true_color);

    %conditional filtering based on some features analysis
    %loading dcube

end