function [x, telap] = amap(y,Phi, sigma, eps)
%AMAP Function to apply the amap algorithm
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps. Returns sparse x and
%   time elapsed
tic
gamma = ones(size(Phi, 2),1);
x = zeros(size(Phi,2),1);
while norm(y-Phi*x) > eps
    s = inv(Phi'*Phi/sigma^2 + inv(diag(gamma)));
    x = s*(Phi'*y)/sigma^2;
    gamma = x.^2;
end
telap = toc;
end

