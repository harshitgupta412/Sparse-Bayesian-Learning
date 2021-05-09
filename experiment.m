rng(4);
% read image
x = double(imread("barbara256.png","png"));

% define patch size and the basis matrix.
patch_size = 8;
U = kron(dctmtx(patch_size)',dctmtx(patch_size)')';

% constants for ista algorithm
Nmax = 100; % number of iterations
lambda = 1; 

counter = zeros(size(x)); %store the count for each pixel, how many time it has been reconstructed
reconstructed_image = zeros(size(x)); %store sum of all the reconstructed patches shifted according to its index

Phi = randn([32, 64]);  % getting the measurement matrix
A = Phi*U;  %The matrix A in ISTA algorithm
alpha = eigs(A'*A,1); %since alpha>eigenvalue of A'A, we can take alpha = max eigenvalue of A'A
        
for i=1:size(x,1)-patch_size+1 
    for j=1:size(x,2)-patch_size+1
        xi = x(i:i+patch_size-1, j:j+patch_size-1); % extract the patch
        yi = Phi * xi(:); %calculate the measurement and vectorising it
        
        % ista
        theta =  zeros(size(xi(:))); %random initialization
        for k = 1:Nmax  
            soft_input = theta + A'*(yi - A*theta)/alpha;
            theta = sign(soft_input) .* max(0, abs(soft_input)- lambda/(2*alpha));
        end
        recovered = (reshape(U*theta, size(xi))); % recover the patch from theta
        
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