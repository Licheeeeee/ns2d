%% 2D incompressible Navier-Stokes solver

%% User settings
% domain size
config.Nx = 40;
config.Ny = 100;
config.dx = 0.01;
config.dy = 0.01;
% operation control
config.Nt = 9000;
config.dt = 1.0;
config.outItvl = 100;
config.useAdv = 1;
% analytical solution
config.v_max = 0.004;
% boundary conditions
config.bcType = [1 1];
config.pBC = [4.8e-4 0];
config.vBC = [0 0];
config.sBC = [1 0];
% initial condition 
config.p0 = 0;
config.u0 = 0;
config.v0 = config.vBC(1);
config.s0 = 0;
% physical parameters
config.nu = 1e-6;
config.kappa = 1e-9;
config.rho = 1000;
config.g = 9.81;
config.Cd = 0.00;
% output
config.makeplot = 0;

%% Bathymetry
config.N = config.Nx * config.Ny;
config.Ng = (config.Nx+2) * (config.Ny+2);
load('domain2D_extended.mat');
% [domain2D] = interpDomain(domain2D,config);
domain = reshape(domain2D,config.N,1);

%% Initialize
[map] = buildMap(domain2D, config);
[data,data2D] = init(map, config);
[data] = enforceBC(data, map, config);

%% Time stepping
Tend = round(config.Nt/config.outItvl);
for tt = 1:config.Nt
    % get ghost velocities
    [data] = updateGhostVel(data, map, config);
    % solution procedure
    [data] = predictor(data, map, config);
    [data] = corrector(data, map, config);
    [data] = updateVelocity(data, map, config);
    [data] = enforceBC(data, map, config);
    [data] = transport(data, map, config);
    % save output
    if round(tt/config.outItvl) == tt/config.outItvl
        [data2D] = saveOut(data, data2D, config, tt/config.outItvl);
    end
%     if nanmax(abs(data.vv)) > 0.99 * config.v_max
%         Tend = floor(tt/config.outItvl);
%         break;
%     end
    % update variables
    data.un = data.uu;
    data.vn = data.vv;
    data.pn = data.pp;
    data.sn = data.ss;
    fprintf('>>>> Time step %d has been executed!\n',tt);
end
% [data2D] = saveOut(data, data2D, config, 1);

