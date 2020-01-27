%% Build map connecting 2D and 1D representation of model domain
function [map] = buildMap(domain2D, config)

    map.cntr = zeros(config.N,1);
    map.iPjc = zeros(config.N,1);
    map.iMjc = zeros(config.N,1);
    map.icjP = zeros(config.N,1);
    map.icjM = zeros(config.N,1);
    map.actv = zeros(config.N,1);
    map.ii = zeros(config.N,1);
    map.jj = zeros(config.N,1);
    map.blockX = zeros(config.Ng,1);
    map.blockY = zeros(config.Ng,1);
    
    kk = 1;
    for jj = 1:config.Ny
        for ii = 1:config.Nx
            if ~isnan(domain2D(ii,jj))
                map.actv(kk) = 1;
            end
            map.ii(kk) = ii;
            map.jj(kk) = jj;
            % center map
            map.cntr(kk) = kk;
            % iM map
            if ii == 1
                map.iMjc(kk) = config.N + 2*config.Nx + config.Ny + jj;
            else
                map.iMjc(kk) = kk - 1;
            end
            % iP map
            if ii == config.Nx
                map.iPjc(kk) = config.N + 2*config.Nx + jj;
            else
                map.iPjc(kk) = kk + 1;
            end
            % jM map
            if jj == 1
                map.icjM(kk) = config.N + config.Nx + ii;
            else
                map.icjM(kk) = kk - config.Nx;
            end
            % jP map
            if jj == config.Ny
                map.icjP(kk) = config.N + ii;
            else
                map.icjP(kk) = kk + config.Nx;
            end
            kk = kk + 1;
        end
    end
    
    % block inactive face
    for kk = 1:config.N
        if map.actv(kk) == 0
            map.blockX(kk) = 1;
            map.blockX(map.iMjc(kk)) = 1;
            map.blockY(kk) = 1;
            map.blockY(map.icjM(kk)) = 1;
        end
    end


end