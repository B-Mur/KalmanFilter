function [xr yr]  = fn_2d(u, ver)

x = u(1); y = u(2);



if ver == 1
    xr = 2*x^2 - 3*x*y + 5*y^2  + sin(x);
    yr = x^3 - y^3 + 5;
elseif ver == 2
    xr = sin(x) * sin(y);
    yr = x^3 - y^2 + cos(x);
elseif ver == 3
    xr = 2 * x;
    yr = 3 * y;
end

    
    
    