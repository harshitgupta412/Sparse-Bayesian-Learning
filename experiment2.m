rng(4);
% read image
x = double(imread("barbara256.png","png"));

%adding noise to get the measurement
y = x + 2*randn(size(x));

% define patch size and the basis matrix.
patch_size = 1;

% constants for ista algorithm
Nmax = 1; % number of iterations
lambda = 1; 
alpha = 1; %since alpha>eigenvalue of A'A, we can take alpha = max eigenvalue of A'A = max eigenvalue of I = 1

counter = zeros(size(x)); %store the count for each pixel, how many time it has been reconstructed
reconstructed_image = zeros(size(x)); %store sum of all the reconstructed patches shifted according to its index
        
for i=1:size(x,1)-patch_size+1 
    for j=1:size(x,2)-patch_size+1
        %since x has sparse representation in 2D dct, we convert the patch of y to dct
        yi = dct2(y(i:i+patch_size-1, j:j+patch_size-1));
        yi = yi(:);
        
        % ista
        theta =  zeros(size(yi)); %random initialization. randn works equally well
        for k = 1:Nmax  
            soft_input = theta + (yi - theta)/alpha;
            theta = sign(soft_input) .* max(0, abs(soft_input)- lambda/(2*alpha));
        end
        recovered = idct2(reshape(theta, [patch_size patch_size])); % recover the patch from theta
        
        % adding patch to image
        reconstructed_image(i:i+patch_size-1, j:j+patch_size-1) = reconstructed_image(i:i+patch_size-1, j:j+patch_size-1) + recovered;
        counter(i:i+patch_size-1, j:j+patch_size-1) = counter(i:i+patch_size-1, j:j+patch_size-1) + 1;
    end
end

% reconstructing and displaying the reconstructed image
x_reconstructed = reconstructed_image ./ counter;
figure, imshow(x/255); title("Original Image");
figure, imshow(x_reconstructed/255); title("Reconstructed Image");

%rmse
rmse = norm(x_reconstructed(:) - x(:))/norm(x(:));
disp("Rmse is "+ string(rmse));