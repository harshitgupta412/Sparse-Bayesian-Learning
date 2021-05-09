% experiment to compare reconstructions for varying sparsity.
% y = x + noise. entries in x is sampled from uniform distribution U(0,1)
% noise variance is 0.1. For each sparsity level, the experiment is
% performed 100 times and averaged.
rng(0);
% sparsity levels
K = [10, 15, 20, 25, 30];
% sparse vector size
total_entry = 100;
% rmse values
rmse_ista = zeros(size(K));
rmse_amap = zeros(size(K));
rmse_sbl = zeros(size(K));
rmse_omp = zeros(size(K));

%time analysis
time_ista = zeros(size(K));
time_amap = zeros(size(K));
time_sbl = zeros(size(K));
time_omp = zeros(size(K));
for index = 1:size(K,2)
    k = K(index);
    %repeat the experiment 100 times for each k
    for re = 1:100
       % constructing the sparse signal x
        numbers = rand(k,1); % non zero values in x
        indices = randperm(total_entry, k); %indices at which the values are non-zero
        x = zeros(total_entry,1);
        x(indices) = numbers;

        % measurements 
        y = x + 0.1*randn(size(x));

        %the matrix A will correspond to Identity
        A = eye(size(x,1));

        % constants for ista algorithm
        Nmax = 100; % number of iterations
        lambda = 1; % regulariser

        % constants for omp, amap, sbl
        sigma = 0.1; % noise standard deviation.
        eps = sigma; % run until ||y-Phi x|| < eps. For gaussian noise, eps = 3*sigma suffices

        % reconstructions
        [tist,telap] = ista(y, A, lambda, Nmax);
        time_ista(index) = time_ista(index) + telap;
        [tsbl, telap] = sbl(y, A,  sigma, eps);
        time_sbl(index) = time_sbl(index) + telap;
        [tomp, telap] = omp(y, A, eps);
        time_omp(index) = time_omp(index) + telap;
        [tamap, telap] = amap(y, A,  sigma, eps);
        time_amap(index) = time_amap(index) + telap;
            
        % rmse
        rmse_omp(index) = rmse_omp(index) + norm(tomp(:) - x(:))/norm(x(:));
        rmse_ista(index) = rmse_ista(index) + norm(tist(:) - x(:))/norm(x(:));
        rmse_amap(index) = rmse_amap(index) + norm(tamap(:) - x(:))/norm(x(:));
        rmse_sbl(index) = rmse_sbl(index) + norm(tsbl(:) - x(:))/norm(x(:));
    end
    rmse_omp(index) = rmse_omp(index)/100;
    rmse_ista(index) = rmse_ista(index)/100;
    rmse_sbl(index) = rmse_sbl(index)/100;
    rmse_amap(index) = rmse_amap(index)/100;
    
    disp("Rmse for k =" + string(k) + " is ");
    disp("OMP: " + string(rmse_omp(index)))
    disp("ISTA: " + string(rmse_ista(index)))
    disp("AMAP: " + string(rmse_amap(index)))
    disp("SBL: " + string(rmse_sbl(index)))
end            

% plotting the rmse values
figure, plot(K/total_entry, rmse_ista, 'g-x', K/total_entry, rmse_sbl, 'b--o', K/total_entry, rmse_amap, 'r:*', K/total_entry, rmse_omp, 'k-.+');
title("Rmse vs Sparsity(number of non zero entries/total entries)")
legend("Ista", "SBL", "AMAP", "OMP", 'Location', 'best');

% plotting the time values
figure, plot(K/total_entry, time_ista/100, 'g-x', K/total_entry, time_sbl/100, 'b--o', K/total_entry, time_amap/100, 'r:*', K/total_entry, time_omp/100, 'k-.+');
title("Time(s) vs Sparsity(number of non zero entries/total entries)")
legend("Ista", "SBL", "AMAP", "OMP", 'Location', 'best');

save('exp2', 'rmse_amap', 'rmse_ista', 'rmse_omp', 'rmse_sbl', 'time_amap', 'time_ista', 'time_omp', 'time_sbl')
