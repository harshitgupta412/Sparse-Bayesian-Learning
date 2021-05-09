function [x,gamma] = amap(y,Phi, sigma, eps)
%AMAP Function to apply the amap algorithm
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps. Returns x and gamma

gamma = ones(1, size(Phi, 2));
x = randn(1, size(Phi,2));
while norm(y-Phi*x) > eps
    s = Phi'*Phi/sigma^2 + inv(diag(gamma));
    x = s\(Phi'*y)/sigma^2;
    gamma = x.^2;
end
end

