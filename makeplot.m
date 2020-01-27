%% Make streamline and scalar figures

r = 1;
fs = 14;
H = 0.2;
y = 1:config.Nx; 
x = 1:config.Ny; 
[X,Y] = meshgrid(x,y);
xstart = []; ystart = [];

s = data2D.s(:,:,Tend);
v = data2D.v(:,:,Tend);
u = data2D.u(:,:,Tend);

aa = map.actv == 0;
s(aa) = -0.1;

for ii = 1:config.Nx
    for jj = 1:4*r:16*r
        xstart = [xstart; ii];
        ystart = [ystart; jj];
    end
    for jj = 17*r:3*r:29*r
        xstart = [xstart; ii];
        ystart = [ystart; jj];
    end
    for jj = 36*r:3*r:48*r
        xstart = [xstart; ii];
        ystart = [ystart; jj];
    end
end

% analytical solution
tVec = [1:Tend] * config.outItvl * config.dt / 3600;
yVec = 0:0.01:H;
dpdx = (config.pBC(1) - config.pBC(2)) / (config.dy * config.Ny);
v_ana = (2*config.rho*config.nu)^(-1) * dpdx .* yVec .* (H - yVec); 

% numerical solution
v_slice = v(11:30,round(0.5*config.Ny));

figure(1);
set(gcf, 'PaperSize', [10,6]);
% streamline
subplot(2,1,1);
streamline(X,Y,v,u,ystart,xstart); hold on; quiver(v,u,'r'); hold off;
xlim([0 size(s,2)*1.1]);
set(gca,'FontSize',fs);
title('Velocity','FontSize',fs);

% scalar
subplot(2,1,2);
imagesc(s);
% contourf(s);
caxis([-0.1 1.0]);
set(gca,'FontSize',fs);
title('Concentration','FontSize',fs);
colorbar;

figure(2);
set(gcf, 'PaperSize', [10,8]);
% velocity profile
subplot(1,2,1);
plot(v_slice,[1:length(v_slice)].*config.dx-0.5*config.dx,'bo');
hold on;
plot(v_ana, yVec,'r-');
grid on;
set(gca, 'YTickLabel', []);
xlabel('Velocity [m/s]','FontSize',fs);
title('Velocity Profile','FontSize',fs);
legend({'Model','Analytical'},'FontSize',fs);

% convergence
subplot(1,2,2);
plot(tVec, squeeze(data2D.v(20,30,1:Tend)));
grid on;
ylabel('Max. Velocity [m/s]','FontSize',fs);
xlabel('Time [h]','FontSize',fs);
title('Velocity convergence','FontSize',fs);






