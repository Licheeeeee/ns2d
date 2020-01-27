%% Function to calculate d(uu)/dx
function [duudx] = ddx(u1, u2, v1, v2, dx)
    duudx = (u1*v1 - u2*v2) / dx;
end