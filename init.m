%% Initialize data arrays
function [data,data2D] = init(map, config)

    data.pp = config.p0 * ones(config.Ng,1);
    data.pn = data.pp;
    
    data.uu = config.u0 * ones(config.Ng,1);
    data.un = data.uu;
    data.us = data.uu;
    data.uy = zeros(config.Ng,1);
    
    data.vv = config.v0 * ones(config.Ng,1);
    data.vn = data.vv;
    data.vs = data.vv;
    data.vx = zeros(config.Ng,1);
    
    data.ugM = NaN * ones(config.Ng,1);
    data.ugP = NaN * ones(config.Ng,1);
    data.vgM = NaN * ones(config.Ng,1);
    data.vgP = NaN * ones(config.Ng,1);
    
    data.ss = config.s0 * ones(config.Ng,1);
%     for kk = 1:config.N
%         if map.actv(kk) == 1
%             data.ss(kk) = config.sBC(1) - map.jj(kk)/config.Ny;
%         end
%     end
    data.sn = data.ss;
    
    data2D.u = zeros(config.Nx,config.Ny,config.Nt/config.outItvl);
    data2D.v = zeros(config.Nx,config.Ny,config.Nt/config.outItvl);
    data2D.p = zeros(config.Nx,config.Ny,config.Nt/config.outItvl);
    data2D.s = zeros(config.Nx,config.Ny,config.Nt/config.outItvl);

end