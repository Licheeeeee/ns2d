%% Convert 1D data into 2D output
function [data2D] = saveOut(data, data2D, config, tt)

    data2D.p(:,:,tt) = reshape(data.pp(1:config.N),config.Nx,config.Ny);
    data2D.u(:,:,tt) = reshape(data.uu(1:config.N),config.Nx,config.Ny);
    data2D.v(:,:,tt) = reshape(data.vv(1:config.N),config.Nx,config.Ny);
    data2D.s(:,:,tt) = reshape(data.ss(1:config.N),config.Nx,config.Ny);
        
end