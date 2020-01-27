%% Enforce boundary conditions
function [data] = enforceBC(data, map, config)

    for kk = 1:config.N
        % blocked face, dp=0, u=0, v=0
        if map.blockX(kk) == 1
            data.uu(kk) = 0.0;
        end
        if map.blockY(kk) == 1
            data.vv(kk) = 0.0;
        end
    end
    
    for kk = 1:config.Nx
        jj = config.Nx*(config.Ny-1) + kk;
        % ym boundary, fixed pressure, zero dv/dy
        if map.actv(kk) == 1
            data.pp(map.icjM(kk)) = config.pBC(1);
            data.uu(map.icjM(kk)) = 0.0;
        else
            data.pp(kk) = 0;
            data.pp(map.icjM(kk)) = 0;
            data.uu(map.icjM(kk)) = 0;
            data.vv(map.icjM(kk)) = 0;
        end
        
        % yp boundary, fixed pressure, velocity from mass conservation
        if map.actv(jj) == 1
            data.pp(map.icjP(jj)) = config.pBC(2);
            data.uu(map.icjP(jj)) = 0.0;
%             data.vv(map.icjP(jj)) = data.vv(map.icjM(kk));
            data.vv(map.icjP(jj)) = data.vv(jj);
        else
            data.pp(jj) = 0;
            data.pp(map.icjP(jj)) = 0;
            data.uu(map.icjP(jj)) = 0;
            data.vv(jj) = 0;
            data.vv(map.icjP(jj)) = 0;
        end
        
    end
    
end

