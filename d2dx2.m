%% Function to calculate d2u/dx2
function [d2udx2] = d2dx2(u, up, um, dx)
    d2udx2 = (up - 2*u + um) / (dx^2);
end