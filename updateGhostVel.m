%% update velocity for ghost cells
function [data] = updateGhostVel(data, map, config)

    for kk = 1:config.N
        if map.blockY(kk) == 1
            data.ugP(kk) = -data.un(kk);
%             data.ugP(map.iMjc(kk)) = -data.un(map.iMjc(kk));
        end
        if map.icjM(kk) <= config.N && map.blockY(map.icjM(kk)) == 1
            data.ugM(kk) = -data.un(kk);
%             data.ugM(map.iMjc(kk)) = -data.un(map.iMjc(kk));
        end
        if map.blockX(kk) == 1
            data.vgP(kk) = -data.vn(kk);
        end
        if map.iMjc(kk) <= config.N && map.blockX(map.iMjc(kk)) == 1
            data.vgM(kk) = -data.vn(kk);
        end
    end


end