function [data] = advection(data, map, config)

    for kk = 1:config.N
        % x advection
        uip = 0.5 * (data.un(kk) + data.un(map.iPjc(kk)));
        uim = 0.5 * (data.un(kk) + data.un(map.iMjc(kk)));
        ujp = 0.5 * (data.un(kk) + data.un(map.icjP(kk)));
        ujm = 0.5 * (data.un(kk) + data.un(map.icjM(kk)));
        
        vjp = 0.5 * (data.vn(kk) + data.vn(map.iPjc(kk)));
        if map.iPjc(kk) <= config.N
            vjm = 0.5 * (data.vn(map.icjM(map.iPjc(kk))) + data.vn(map.icjM(kk)));
        else
            vjm = 0.5 * data.vn(map.icjM(kk));
        end

        if map.blockX(map.icjM(kk)) == 1
            ujm = 0;
            vjm = 0;
        end
        if map.blockX(map.icjP(kk)) == 1
            ujp = 0;
            vjp = 0;
        end
        
        data.us(kk) = data.us(kk) - (config.dt/config.dx) * (uip*uip - uim*uim) - ...
            (config.dt/config.dy) * (ujp*vjp - ujm*vjm);
        
        % y advection
        vip = 0.5 * (data.vn(kk) + data.vn(map.iPjc(kk)));
        vim = 0.5 * (data.vn(kk) + data.vn(map.iMjc(kk)));
        vjp = 0.5 * (data.vn(kk) + data.vn(map.icjP(kk)));
        vjm = 0.5 * (data.vn(kk) + data.vn(map.icjM(kk)));
        
        uip = 0.5 * (data.un(kk) + data.un(map.icjP(kk)));
        if map.iMjc(kk) <= config.N
            uim = 0.5 * (data.un(map.icjP(map.iMjc(kk))) + data.un(map.iMjc(kk)));
        else
            uim = 0.5 * data.un(map.iMjc(kk));
        end
        
        if map.blockY(map.iMjc(kk)) == 1
            uim = 0;
            vim = 0;
        end
        if map.blockY(map.iPjc(kk)) == 1
            uip = 0;
            vip = 0;
        end
        
        data.vs(kk) = data.vs(kk) - (config.dt/config.dx) * (vip*uip - vim*uim) - ...
            (config.dt/config.dy) * (vjp*vjp - vjm*vjm);
        
    end
    
    % ghost cells
    for kk = 1:config.Nx
        mm = map.icjM(kk);
        jj = config.Nx*(config.Ny-1) + kk;
        if map.blockX(kk) == 1
            vip = 0;
        else
            vip = 0.5 * (data.vn(mm) + data.vn(mm+1));
        end
        if map.blockX(map.iMjc(kk)) == 1
            vim = 0;
        else            
            vim = 0.5 * (data.vn(mm) + data.vn(mm-1));
        end
        uip = 0.5 * data.un(kk);
        uim = 0.5 * data.un(map.iMjc(kk));
        
        
        vjp = 0.5 * (data.vn(map.icjM(kk)) + data.vn(kk));
%         vjm = 0.5 * (data.vn(map.icjM(kk)) + data.vn(jj));
        vjm = 0.5 * (data.vn(map.icjM(kk)) + data.vn(map.icjM(kk)));
        data.vs(map.icjM(kk)) = data.vs(map.icjM(kk)) - ...
            (config.dt/config.dx) * (vip*uip - vim*uim) - ...
            (config.dt/config.dy) * (vjp*vjp - vjm*vjm);
        data.us(map.icjM(kk)) = 0.0;
    end


end