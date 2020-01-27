%% Solve the Poisson equation to correct pressure
function [data] = corrector(data, map, config)
    Na = sum(map.actv(1:config.N));
    A = sparse(Na,Na);
    B = zeros(Na,1);
    row = 0;
    for kk = 1:config.N
        if map.actv(kk) == 1
            row = row + 1;
            % build linear system
            A(row,row) = 0.0;
            B(row) = (config.rho*config.dt/config.dx)*(data.us(kk)-data.us(map.iMjc(kk))) + ...
                (config.rho*config.dt/config.dy)*(data.vs(kk)-data.vs(map.icjM(kk)));
            if map.blockX(kk) == 0 && map.iPjc(kk) <= config.N
                A(row,row) = A(row,row) - 1/config.dx^2;
                A(row,row+1) = 1/config.dx^2;
            end
            if map.blockX(map.iMjc(kk)) == 0 && map.iMjc(kk) < config.N
                A(row,row) = A(row,row) - 1/config.dx^2;
                A(row,row-1) = 1/config.dx^2;
            end
            if map.blockY(kk) == 0
                A(row,row) = A(row,row) - 1/config.dy^2;
                if map.icjP(kk) <= config.N
                    dist = sum(map.actv(kk:kk+config.Nx)) - 1;
                    A(row,row+dist) = 1/config.dy^2;
                else
                    B(row) = B(row) - config.pBC(2)/config.dy^2;
                end
            end
            if map.blockY(map.icjM(kk)) == 0
                A(row,row) = A(row,row) - 1/config.dy^2;
                if map.icjM(kk) <= config.N
                    dist = sum(map.actv(kk-config.Nx:kk)) - 1;
                    A(row,row-dist) = 1/config.dy^2;
                else
                    B(row) = B(row) - config.pBC(1)/config.dy^2;
                end
            end
        end
    end
    
    % solve the linear system
    x = A \ B;
    
    % save the matrix for debugging
    data.A = A;
    data.B = B;
    data.x = x;
    
    row = 0;
    for kk = 1:config.N
        if map.actv(kk) == 1
            row = row + 1;
            data.pp(kk) = x(row);
        else
            data.pp(kk) = config.p0;
        end
    end
    
    for kk = 1:config.Nx
        jj = config.Nx*(config.Ny-1)+kk;
        if map.actv(kk) == 1
            data.pp(map.icjM(kk)) = config.pBC(1);
        end
        if map.actv(jj) == 1
            data.pp(map.icjP(jj)) = config.pBC(2);
        end
    end
    
    
    

end