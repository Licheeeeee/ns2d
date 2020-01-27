%% Predictor step to get us, vs
function [data] = predictor(data, map, config)
    
    % calculate us, vs
    for kk = 1:config.N
        data.us(kk) = data.un(kk);
        data.vs(kk) = data.vn(kk);
    end
    for kk = 1:config.Nx
        data.us(map.icjM(kk)) = 0;
        data.vs(map.icjM(kk)) = data.vn(map.icjM(kk));
    end        
    
    % diffusion
    [data] = diffusion(data, map, config);
    
    % advection
    if config.useAdv == 1
        [data] = advection(data, map, config);
    end
    
    % friction
    for kk = 1:config.N
        vel = sqrt(data.us(kk)^2 + data.vs(kk)^2);
        if vel > 0
            data.us(kk) = data.us(kk) / (1 + 0.5*config.dt*config.Cd*vel);
            data.vs(kk) = data.vs(kk) / (1 + 0.5*config.dt*config.Cd*vel);
        end
    end

end