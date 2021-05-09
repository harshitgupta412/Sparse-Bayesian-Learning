rng(4);
% values of number of measurements
M = [8, 16, 32, 48, 64];
% rmse values
rmse_ista = zeros(size(M));
rmse_amap = zeros(size(M));
rmse_sbl = zeros(size(M));
rmse_omp = zeros(size(M));

%time analysis
time_ista = zeros(size(M));
time_amap = zeros(size(M));
time_sbl = zeros(size(M));
time_omp = zeros(size(M));

% read image
x = double(imread("barbara.png","png"));

% constants for ista algorithm
Nmax = 100; % number of iterations
lambda = 1; % regulariser

% constants for omp, amap, sbl
eps = 0.01; % run until ||y-Phi x|| < eps
sigma = 0.01; % noise standard deviation

for index = 1:size(M,2)
    % Number of measurements
    m = M(index);
    % define patch size and the basis matrix.
    patch_size = 8;
    % dct basis
    U = kron(dctmtx(patch_size)',dctmtx(patch_size)')';
    Phi = randn([m, 64]);  % getting the measurement matrix
    A = Phi*U;  %The matrix Phi in algorithms. y = \Phi x

    counter = zeros(size(x)); %store the count for each pixel, how many time it has been reconstructed
    ista_cons = zeros(size(x)); %store sum of all the reconstructed patches from ista shifted according to its index 
    sbl_cons = zeros(size(x)); %store sum of all the reconstructed patches from sbl shifted according to its index
    amap_cons = zeros(size(x)); %store sum of all the reconstructed patches from amap shifted according to its index
    omp_cons = zeros(size(x)); %store sum of all the reconstructed patches from omp shifted according to its index

    for i=1:size(x,1)-patch_size+1 
        disp(i);
        for j=1:size(x,2)-patch_size+1
            xi = x(i:i+patch_size-1, j:j+patch_size-1); % extract the patch
            yi = Phi * xi(:); %calculate the measurement and vectorising it

            % reconstructions
            [tist,telap] = ista(yi, A, lambda, Nmax);
            time_ista(index) = time_ista(index) + telap;
            
            [tsbl, telap] = sbl(yi, A,  sigma, eps,0);
            time_sbl(index) = time_sbl(index) + telap;
            
            
            if m < 64 % for m=64, the matrix is not invertible in amap. For omp it takes a lot of time
               [tomp, telap] = omp(yi, A, eps,0);
                time_omp(index) = time_omp(index) + telap;
                [tamap, telap] = amap(yi, A,  sigma, eps,0);
                time_amap(index) = time_amap(index) + telap;
            end
            
            % adding patch to reconstructed image
            ista_cons(i:i+patch_size-1, j:j+patch_size-1) = reshape(U*tist, size(xi)) + ista_cons(i:i+patch_size-1, j:j+patch_size-1);
            sbl_cons(i:i+patch_size-1, j:j+patch_size-1) = reshape(U*tsbl, size(xi)) + sbl_cons(i:i+patch_size-1, j:j+patch_size-1);
            amap_cons(i:i+patch_size-1, j:j+patch_size-1) = reshape(U*tamap, size(xi)) + amap_cons(i:i+patch_size-1, j:j+patch_size-1);
            omp_cons(i:i+patch_size-1, j:j+patch_size-1) = reshape(U*tomp, size(xi)) + omp_cons(i:i+patch_size-1, j:j+patch_size-1);
            counter(i:i+patch_size-1, j:j+patch_size-1) = counter(i:i+patch_size-1, j:j+patch_size-1) + 1;
        end
    end

    % reconstructing  the image by dividing by the count
    ista_cons = ista_cons ./ counter;
    sbl_cons = sbl_cons ./ counter;
    amap_cons = amap_cons ./ counter;
    omp_cons = omp_cons ./ counter;
        
    % displaying and saving the reconstructions
    figure;
    subplot(2,2,1); imshow(omp_cons/255); title("OMP Reconstruction");
    subplot(2,2,2); imshow(ista_cons/255); title("ISTA Reconstruction");
    subplot(2,2,3); imshow(sbl_cons/255); title("SBL Reconstruction");
    subplot(2,2,4); imshow(amap_cons/255); title("AMAP Reconstruction");
    saveas(gcf, string(m) + '.png');

    %rmse values
    rmse_omp(index) = norm(omp_cons(:) - x(:))/norm(x(:));
    rmse_ista(index) = norm(ista_cons(:) - x(:))/norm(x(:));
    rmse_amap(index) = norm(amap_cons(:) - x(:))/norm(x(:));
    rmse_sbl(index) = norm(sbl_cons(:) - x(:))/norm(x(:));
    disp("Rmse for m =" + string(m) + " is ");
    disp("OMP: " + string(rmse_omp(index)))
    disp("ISTA: " + string(rmse_ista(index)))
    disp("AMAP: " + string(rmse_amap(index)))
    disp("SBL: " + string(rmse_sbl(index)))
end

% calculating average time elapsed per patch = Total time/number of patches
time_ista = time_ista / (size(x,1)-patch_size+1)^2;
time_omp = time_omp / (size(x,1)-patch_size+1)^2;
time_amap = time_amap / (size(x,1)-patch_size+1)^2;
time_sbl = time_sbl / (size(x,1)-patch_size+1)^2;
    
% plotting the rmse values
figure, plot(M, rmse_ista, 'g-x', M, rmse_sbl, 'b--o', M, rmse_amap, 'r:*', M, rmse_omp, 'k-.+');
title("Rmse vs M")
legend("Ista", "SBL", "AMAP", "OMP", 'Location', 'best');

save('barbara', 'rmse_amap', 'rmse_ista', 'rmse_omp', 'rmse_sbl', 'time_amap', 'time_ista', 'time_omp', 'time_sbl')
