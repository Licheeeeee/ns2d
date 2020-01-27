%% Scalar transport
function [data] = transport(data, map, config)

    for kk = 1:config.N
        if map.actv(kk) == 1
            % advective transport
            Qxp = data.uu(kk) * config.dy;
            Qxm = data.uu(map.iMjc(kk)) * config.dy;
            Qyp = data.vv(kk) * config.dx;
            Qym = data.vv(map.icjM(kk)) * config.dy;
            V = config.dx * config.dy;
            
            data.ss(kk) = data.sn(kk) * V;
            if Qxp > 0
                data.ss(kk) = data.ss(kk) - config.dt * Qxp * data.sn(kk);
            elseif Qxp < 0
                data.ss(kk) = data.ss(kk) - config.dt * Qxp * data.sn(map.iPjc(kk));
            end
            if Qxm > 0
                data.ss(kk) = data.ss(kk) + config.dt * Qxm * data.sn(map.iMjc(kk));
            elseif Qxm < 0
                data.ss(kk) = data.ss(kk) + config.dt * Qxm * data.sn(kk);
            end
            if Qyp > 0
                data.ss(kk) = data.ss(kk) - config.dt * Qyp * data.sn(kk);
            elseif Qyp < 0
                data.ss(kk) = data.ss(kk) - config.dt * Qyp * data.sn(map.icjP(kk));
            end
            if Qym > 0
                data.ss(kk) = data.ss(kk) + config.dt * Qym * data.sn(map.icjM(kk));
            elseif Qym < 0
                data.ss(kk) = data.ss(kk) + config.dt * Qym * data.sn(kk);
            end
            
            % diffusive transport 
            if map.blockX(kk) == 0
                data.ss(kk) = data.ss(kk) + (config.dt * config.kappa * config.dx) * ...
                    (data.sn(map.iPjc(kk)) - data.sn(kk)) / config.dy;
            end
            if map.iMjc(kk) <= config.N && map.blockX(map.iMjc(kk)) == 0
                data.ss(kk) = data.ss(kk) + (config.dt * config.kappa * config.dx) * ...
                    (data.sn(kk) - data.sn(map.iMjc(kk))) / config.dy;
            end
            if map.blockY(kk) == 0
                data.ss(kk) = data.ss(kk) + (config.dt * config.kappa * config.dy) * ...
                    (data.sn(map.icjP(kk)) - data.sn(kk)) / config.dx;
            end
            if map.icjM(kk) <= config.N && map.blockY(map.icjM(kk)) == 0
                data.ss(kk) = data.ss(kk) + (config.dt * config.kappa * config.dy) * ...
                    (data.sn(kk) - data.sn(map.icjM(kk))) / config.dx;
            elseif map.icjM(kk) > config.N && map.actv(kk) == 1
                data.ss(kk) = data.ss(kk) + (config.dt * config.kappa * config.dy) * ...
                    (data.sn(kk) - data.sn(map.icjM(kk))) / config.dx;
            end
            
            % divided by volume to get concentration
            data.ss(kk) = data.ss(kk) / V;
            
            % scalar limiter to remove local maxima/minima
            if map.iPjc(kk) <= config.N && map.iMjc(kk) <= config.N && ...
                    map.icjP(kk) <= config.N && map.icjM(kk) <= config.N
                actv4 = [map.actv(map.iPjc(kk)) map.actv(map.iMjc(kk)) map.actv(map.icjP(kk)) ...
                    map.actv(map.icjM(kk))];
                % boundary scalar
                if map.actv(kk) == 1 && map.icjM(kk) > config.N
                    actv4(4) = 1;
                end
                S4 = [data.sn(map.iPjc(kk)) data.sn(map.iMjc(kk)) data.sn(map.icjP(kk)) ...
                    data.sn(map.icjM(kk))];
                aa = actv4 == 1;
                S4 = S4(aa);
                if ~isempty(S4)
                    if data.ss(kk) > max(S4)
                        data.ss(kk) = max(S4);
                    end
                    if data.ss(kk) < min(S4)
                        data.ss(kk) = min(S4);
                    end
                end
            end
            
        else
            data.ss(kk) = 0.0;
        end
    end
    
    for kk = 1:config.Nx
        data.ss(map.icjM(kk)) = config.sBC(1);
    end
        
    


end