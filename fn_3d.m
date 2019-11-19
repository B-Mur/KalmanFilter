function [xr, yr, zr]  = fn_3d(u, ver)

x = u(1); y = u(2); z = u(3);

if ver == 1
    xr = 2*x^3 + (3*x*y)^2 + 5*y^2 - cos(y);
    yr = x^3 - y^3 + 5;
    zr = y*z + abs(x)^(1/3);
elseif ver == 2
    xr = y^3 + (x*z)^2 +  sin(y);
    yr = z^3 - x^3 + cos(x);
    zr = y*z + (x)*cos(y);   
    
elseif ver == 3
    xr = 2*y;
    yr = 2*x;
    zr = z*5;
end
