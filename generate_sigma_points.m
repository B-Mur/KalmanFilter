function [sig_pts]  = generate_sigma_points(x, P)

d = size(P, 1);
N = 2*d; 
X = randn(d, d);
X = [X; -X]; 

C = cov(X) * (N-1)/N;
X = C^(-1/2) * X';
X = sqrtm(P) * X;

sig_pts = X + x';

if cov(sig_pts')*(N-1)/N - P>2*eps
    disp('sigma points covariance does not match P');
    disp('if you made it here, you have some serious problems');
    keyboard; 
elseif mean(sig_pts')-x >2*eps
    disp('sigma points mean does not match x');
    disp('if you made it here, you have some serious problems');

    keyboard; 
end

