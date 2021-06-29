function [x, telap] = omp(y,Phi,eps, Nmax)
%OMP Function to apply the sbl (EM) algorithm
%   y: measurements, Phi: measurement matrix, The
%   algorithm is executed till ||residual|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
tic
r = y; % residual
Ts = []; % support
norm_Phi = normc(Phi);

if Nmax == 0
    while norm(r) > eps
        dot_products= abs(r'*norm_Phi); 
        [~, argmax] = max(dot_products); % this gives us the
        % maximum dot product of r with the columns
        Ts = [Ts argmax]; % adding to support

        subA = Phi(:, Ts); % choosing the submatrix corresponding to support
        theta = pinv(subA)*y; % projected distance minimizer theta
        r = y - subA*theta; % remaining residual
    end
else
    Nmax = min(Nmax, size(y,1));
    for re = 1:Nmax
        dot_products= abs(r'*norm_Phi); 
        [~, argmax] = max(dot_products); % this gives us the
        % maximum dot product of r with the columns
        Ts = [Ts argmax]; % adding to support

        subA = Phi(:, Ts); % choosing the submatrix corresponding to support
        theta = pinv(subA)*y; % projected distance minimizer theta
        r = y - subA*theta; % remaining residual
    end
end

x = zeros(size(Phi, 2), 1); 
x(Ts) = theta; % this is the sparse vector we obtained 
telap = toc;
end
