%% Update velocity from pressure
function [data] = updateVelocity(data, map, config)
    
    for kk = 1:config.N
        if map.blockX(kk) == 0
            data.uu(kk) = data.us(kk) - (config.dt/(config.dx*config.rho)) * (data.pp(map.iPjc(kk)) - data.pp(kk));
        else
            data.uu(kk) = 0.0;
        end
        if map.blockY(kk) == 0
            data.vv(kk) = data.vs(kk) - (config.dt/(config.dy*config.rho)) * (data.pp(map.icjP(kk)) - data.pp(kk));
        else
            data.vv(kk) = 0.0;
        end
        if abs(data.uu(kk)) < 1e-10
            data.uu(kk) = 0.0;
        end
        if abs(data.vv(kk)) < 1e-10
            data.vv(kk) = 0.0;
        end
    end
        
    for kk = 1:config.Nx
        mm = map.icjM(kk);
        data.uu(mm) = 0.0;
        if map.actv(kk) == 1
            data.vv(mm) = data.vs(mm) - (config.dt/(config.dy*config.rho)) * (data.pp(kk) - data.pp(mm));
        else
            data.vv(mm) = 0.0;
        end
    end

end