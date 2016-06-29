% =========================================================================
% Purpose - Numerical solutions to 2D model of Marshall and Radko 2003 
%           Residual mean solution for the ACC and its assosciated
%           overturning circulation
% Author - Dhruv Balwada
% Date   - 2015
% =========================================================================
% Note   - Solved in finite difference form and not using the method of
% characteristics that is used in the paper.
% =========================================================================

clear all 
close all 

dz = 5; 
dy = 50*1000; 

% 10^-4 is chosen based on manual entry, ...
% proper way would be to calculate CFL in the compute loop

CFL = dz/dy/10^(-4);  

if CFL >1
    disp('CFL too high, reduce grid size.')
    exit
end

z = -[hm:dz:2500]; 
y = [0:dy:2000*1000]; 

psi_bar  = nan(length(y),length(z));
psi_star = nan(length(y),length(z));
psi_res  = nan(length(y),length(z));
b        = nan(length(y),length(z)); 
wc       = b;

% Parameters and forcings
hm    = 200;         % Mixed layer depth
tau_o = 2*10^(-4);   % surface wind stress scale
Ly    = 2000*1000;   % Meridional scale
Lx    = 21000*1000;  % Zonal Width 
tau   = tau_o*(0.3 + sin(pi*y/Ly)); 
delbo = 7*10^(-3);   % buoyancy change across ACC 
bo    = delbo*y/Ly;  
Bo    = 2*10^(-9);   % Net byouancy forcing 
B     = -Bo*sin(2*pi*y/Ly); 
f     = -10^(-4); 
beta  = 10^(-11); 
k     = 10^6;        % eddy parameter

id_hm = find(z==-hm);
psi_res(:,id_hm) = B/(delbo/Ly); % choose this because the gradient is a constant here
                                 % otherwise take the proper db/dy
                                 % gradient; 

psi_res_b = B/(delbo/Ly);        % The residual circulation based on the
                                 % mixed layer buoyancy budget

b(:,1) = bo;    % specify the surface buoyancy
b(1,:) = bo(1); % southern wall has same bo as at the surface over there

for j=2:length(z)-1 % loop from surface to depth 
        idb = find(bo<=b(1,j-1),1,'last'); % find the buoyancy at the base of ML
        if ~isempty(idb)
            wc(1,j-1) = -sqrt(-tau(1)/f/k - psi_res_b(idb)/k);     
            b(1,j) = b(1,j-1) + (-dz)/wc(1,j-1)*(b(2,j-1)-b(1,j-1))/dy;
        else
            b(1,j)= 0;
        end
    for i =2:length(y)-1 % loop over the meridional extent
        idb = find(bo<=b(i,j-1),1,'last'); % find the buoyancy at the at the base of ML
        if ~isempty(idb)
            wc(i,j-1) = -sqrt(-tau(i)/f/k - psi_res_b(idb)/k);     
            b(i,j) = b(i,j-1) - (-dz)/wc(i,j-1)*(b(i+1,j-1)-b(i-1,j-1))/2/dy;
        else
            b(i,j)=0;
        end
    end
        idb = find(bo<=b(end,j-1),1,'last'); % find the buoyancy at the at the base of ML
        if ~isempty(idb)
            wc(end,j-1) = -sqrt(-tau(end)/f/k - psi_res_b(idb)/k);     
            b(end,j) = b(end,j-1) - (-dz)/wc(end,j-1)*(b(end,j-1)-b(end-1,j-1))/dy;
        else
            b(end,j) = 0;
        end
    
end

%% Plot it. 

h= contourf(y,z,b',20);
caxis([min(bo) max(bo)])
colorbar