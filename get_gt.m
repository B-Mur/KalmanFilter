function [t_norm_dist2d, t_norm_dist3d, t_uni_dist2d, t_uni_dist3d]  = get_gt(norm_dist2, uni_dist2, norm_dist3, uni_dist3, ver)

N = size(norm_dist2, 1)
t_norm_dist2d = zeros(N, 2); t_norm_dist3d = zeros(N, 3);
t_uni_dist2d  = zeros(N, 2); t_uni_dist3d  = zeros(N, 3);

for i = 1:size(norm_dist2, 1)
  u = norm_dist2(i,:); 
  [x2, y2] = fn_2d(u, ver);
  u = norm_dist3(i,:); 
  [x3, y3, z3] = fn_3d(u, ver);
    
  % Normal Distributions
  t_norm_dist2d(i,:) = [x2, y2];
  t_norm_dist3d(i,:) = [x3, y3, z3];

  % Uniform Distribution
  u = uni_dist2(i,:);  
  [x2, y2] = fn_2d(u(1:2), ver);
  u = uni_dist3(i,:); 
  [x3, y3, z3] = fn_3d(u, ver);
  
  t_uni_dist2d(i,:) = [x2, y2];
  t_uni_dist3d(i,:) = [x3, y3, z3];
  if (mod(i,100000) == 0)
    fprintf('Iteration : %i \n', i)
  end

end
