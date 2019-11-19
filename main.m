
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%

N = 100000;
ver = 1; % WHICH FUNCIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
syms f2(x,y) g2(x,y) f3(x,y,z) g3(x,y,z) h3(x,y,z)

if ver == 1
f2(x,y) = 2*x^2 - 3*x*y + 5*y^2 + sin(x);
g2(x,y) = x^3 - y^3 + 5;

f3(x,y,z) = 2*x^3 + (3*x*y)^2 + 5*y^2 - cos(y);
g3(x,y,z) = x^3 - y^3 + 5;
h3(x,y,z) = y*z + abs(x)^(1/3);

elseif ver == 2
f2(x,y) = sin(x) * sin(y);
g2(x,y) = x^3 - y^2 + cos(x);

f3(x,y,z) = y^3 + (x*z)^2 +  sin(y);
g3(x,y,z) = z^3 - x^3 + cos(x);
h3(x,y,z) = y*z + (x)*cos(y);

elseif ver == 3
    
f2(x,y) = 2 * x;
g2(x,y) = 3 * y;
    
f3(x,y,z) = 2*x;
g3(x,y,z) = 3*y;
h3(x,y,z) = 5*z;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate distributions

norm_dist2 = [mvnrnd(.1,.2,N), mvnrnd(.6,2,N)] ; % Gaussian - Distribution 1
uni_dist2 =  [unifrnd(2, 3, [1, N]); unifrnd(3, 7, [1, N])]' ;% Uniform - Distribution 2

norm_dist3 = [mvnrnd(.05,.05,N), mvnrnd(.85, 2, N), mvnrnd(.9,6,N)] ; % Gaussian - Distribution 1
uni_dist3 =  [unifrnd(2, 4, [1, N]); unifrnd(3, 4, [1, N]); unifrnd(3, 3.5, [1, N])]' ;   % Uniform - Distribution 2


mu_norm2 = mean(norm_dist2);
mu_uni2  = mean(uni_dist2);
cov_norm2 = cov(norm_dist2) * (N-1)/N;
cov_uni2 = cov(uni_dist2) * (N-1)/N;

mu_norm3 = mean(norm_dist3);
mu_uni3  = mean(uni_dist3);
cov_norm3 = cov(norm_dist3) * (N-1)/N;
cov_uni3 = cov(uni_dist3) * (N-1)/N;
 


%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 
%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform
    % initialize


[t_norm_dist2d, t_norm_dist3d, t_uni_dist2d, t_uni_dist3d]  = get_gt(norm_dist2, uni_dist2, norm_dist3, uni_dist3, ver);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the ground truth means and covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gt_mean2d_norm = mean(t_norm_dist2d);
gt_cov2d_norm = cov(t_norm_dist2d) * ((N-1)/N);

gt_mean3d_norm = mean(t_norm_dist3d);
gt_cov3d_norm = cov(t_norm_dist3d) * ((N-1)/N);

gt_mean2d_uni = mean(t_uni_dist2d);
gt_cov2d_uni = cov(t_uni_dist2d) * ((N-1)/N);

gt_mean3d_uni = mean(t_uni_dist3d);
gt_cov3d_uni = cov(t_uni_dist3d) * ((N-1)/N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Jacobians
jacob_2d = jacobian([f2, g2], [x, y]);
jacob_3d = jacobian([f3,g3, h3], [x, y, z]);

%%%%%%%%%%%%%%%%%%
%     2D         %
%%%%%%%%%%%%%%%%%%
mu_norm2d = mean(norm_dist2);
Si_norm = cov(norm_dist2) * ((N-1)/N); 
grad_norm = double(jacob_2d(mu_norm2d(1), mu_norm2d(2)));
[new_x, new_y] = fn_2d(mu_norm2d, ver);
new_mu_2norm = [new_x; new_y];
new_Si_2norm = grad_norm * Si_norm * grad_norm';

norm_2d_mahal_lin = bryce_mahal(new_mu_2norm, new_Si_2norm, gt_mean2d_norm', gt_cov2d_norm)

mu_unif2d = mean(uni_dist2);
Si_unif =  cov(uni_dist2) * ((N-1)/N);
grad_uni = double(jacob_2d(mu_unif2d(1), mu_unif2d(2)));
[new_x, new_y] = fn_2d(mu_unif2d, ver);
new_mu_2uni = [new_x; new_y];
new_Si_2uni = grad_uni * Si_unif * grad_uni';

uni_2d_mahal_lin = bryce_mahal(new_mu_2uni, new_Si_2uni, gt_mean2d_uni', gt_cov2d_uni)


%%%%%%%%%%%%%%%%%%
%     3D         %
%%%%%%%%%%%%%%%%%%

mu_norm3d = mean(norm_dist3);
Si_norm = cov(norm_dist3) * ((N-1)/N); 
grad_norm = double(jacob_3d(mu_norm3d(1), mu_norm3d(2), mu_norm3d(3)));
[x, y, z] = fn_3d(mu_norm3d, ver); % Need to output 3 things, not one, dummy
new_mu_3norm = [x, y, z];
new_Si_3norm = grad_norm * Si_norm * grad_norm';

norm_3d_mahal_lin = bryce_mahal(new_mu_3norm', new_Si_3norm, gt_mean3d_norm', gt_cov3d_norm)


mu_unif3 = mean(uni_dist3);
Si_unif =  cov(uni_dist3) * ((N-1)/N);
grad_uni = double(jacob_3d(mu_unif3(1), mu_unif3(2), mu_unif3(3)));
[x, y, z] = fn_3d(mu_unif3, ver); % Need to output 3 things, not one, dummy
new_mu_3uni = [x, y, z];
new_Si_3uni = grad_uni * Si_unif * grad_uni';

uni_3d_mahal_lin = bryce_mahal(new_mu_3uni', new_Si_3uni, gt_mean3d_uni', gt_cov3d_uni)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unscented Kalman Filter : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate sigmas points - norm
sig_pts_norm_2d = generate_sigma_points(mu_norm2, cov_norm2);

sig_pts_norm_3d = generate_sigma_points(mu_norm3, cov_norm3);

% generate sigmas points - uni
sig_pts_uni_2d = generate_sigma_points(mu_uni2, cov_uni2);

sig_pts_uni_3d = generate_sigma_points(mu_uni3, cov_uni3);

% apply function to each of the points

for i = 1:size(sig_pts_norm_2d, 2)
    
   [x, y] = fn_2d(sig_pts_norm_2d(:,i), ver);
   sig_norm2d_res(i,:) = [x, y];
   
   [x, y] = fn_2d(sig_pts_uni_2d(:,i), ver);
   sig_uni2d_res(i,:) =[x, y] ;

  
   
end

for i = 1:size(sig_pts_norm_3d, 2)
       [x, y, z] = fn_3d(sig_pts_norm_3d(:,i), ver);
   sig_norm3d_res(i,:) = [x, y, z];
   
   [x, y, z] = fn_3d(sig_pts_uni_3d(:,i), ver);
   sig_uni3d_res(i,:) = [x, y, z];
   
end

NN2 = size(sig_pts_norm_2d, 2);
NN3 = size(sig_pts_norm_3d, 2);

mean_unsc_norm2d = mean(sig_norm2d_res);
sigm_unsc_norm2d = cov(sig_norm2d_res) * (NN2-1)/NN2;
norm_2d_mahal_unsc = bryce_mahal(mean_unsc_norm2d', sigm_unsc_norm2d, gt_mean2d_norm', gt_cov2d_norm)


mean_unsc_uni2d = mean(sig_uni2d_res);
sigm_unsc_uni2d = cov(sig_uni2d_res) * (NN2-1)/NN2;
uni_2d_mahal_unsc = bryce_mahal(mean_unsc_uni2d', sigm_unsc_uni2d, gt_mean2d_uni', gt_cov2d_uni)

mean_unsc_norm3d = mean(sig_norm3d_res);
sigm_unsc_norm3d = cov(sig_norm3d_res) * (NN3-1)/NN3;
norm_3d_mahal_unsc = bryce_mahal(mean_unsc_norm3d', sigm_unsc_norm3d, gt_mean3d_norm', gt_cov3d_norm)


mean_unsc_uni3d = mean(sig_uni3d_res);
sigm_unsc_uni3d = cov(sig_uni3d_res) * (NN3-1)/NN3;
uni_3d_mahal_unsc = bryce_mahal(mean_unsc_uni3d', sigm_unsc_uni3d, gt_mean3d_uni', gt_cov3d_uni)







