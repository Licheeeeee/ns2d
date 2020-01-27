%% Compute diffusion
function [data] = diffusion(data, map, config)

    for kk = 1:config.N
        % x-diffusion
        data.us(kk) = data.us(kk) + config.dt * config.nu * ...
                d2dx2(data.un(kk),data.un(map.iPjc(kk)),data.un(map.iMjc(kk)),config.dx);
        if map.blockX(kk) == 0
            if map.icjP(kk) <= config.N && map.blockX(map.icjP(kk)) == 1
                ujP = -data.un(kk);
            else
                ujP = data.un(map.icjP(kk));
            end
            if map.icjM(kk) <= config.N && map.blockX(map.icjM(kk)) == 1
                ujM = -data.un(kk);
            else
                ujM = data.un(map.icjM(kk));
            end
            data.us(kk) = data.us(kk) + config.dt * config.nu * ...
                d2dx2(data.un(kk),ujP,ujM,config.dy);
        end
        % y-diffusion
        data.vs(kk) = data.vs(kk) + config.dt * config.nu * ...
                d2dx2(data.vn(kk),data.vn(map.icjP(kk)),data.vn(map.icjM(kk)),config.dy);
        if map.blockY(kk) == 0
            if map.iPjc(kk) <= config.N && map.blockY(map.iPjc(kk)) == 1
                viP = -data.vn(kk);
            else
                viP = data.vn(map.iPjc(kk));
            end
            if map.iMjc(kk) <= config.N && map.blockY(map.iMjc(kk)) == 1
                viM = -data.vn(kk);
            else
                viM = data.vn(map.iMjc(kk));
            end
            data.vs(kk) = data.vs(kk) + config.dt * config.nu * ...
                d2dx2(data.vn(kk),viP,viM,config.dx);
        end
        
    end
    
    % ghost cells
    for kk = 1:config.Nx
        mm = map.icjM(kk);
        jj = config.Nx*(config.Ny-1)+kk;
        if map.blockX(kk) == 1
            viP = -data.vn(mm);
        else
            viP = data.vn(mm+1);
        end
        if map.blockX(map.iMjc(kk)) == 1
            viM = -data.vn(mm);
        else            
            viM = data.vn(mm-1);
        end
        % periodic BC
%         data.vs(mm) = data.vs(mm) + config.dt * config.nu * ...
%                 (d2dx2(data.vn(mm),viP,viM,config.dx) + ...
%                 d2dx2(data.vn(mm),data.vn(kk),data.vn(jj),config.dy));
        data.vs(mm) = data.vs(mm) + config.dt * config.nu * ...
                (d2dx2(data.vn(mm),viP,viM,config.dx) + ...
                d2dx2(data.vn(mm),data.vn(kk),data.vn(mm),config.dy));
    end


end