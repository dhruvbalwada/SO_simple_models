% This script reproduces the results of the paper Badin and Williams 2010
% A model which incorporates from Mixed layer physics to compare its
% interactions with the overturning circulation. Slight improvement over
% the Marshall and Radko 2003 work, which assumed a buoyancy relaxation
% condition at the surface.

% parameters
close all 
clear all 

global f Lx Ly H hm cpa rhoa rhow cd a b c eps tau0 K k sbc Ah y z tau Tair uw Pa

f = -10^(-4);
Lx = 21*10^6;
Ly = 2*10^6;
H = 1000;
hm = 200;          % mixed layer depth
cpa = 1008;
rhoa = 1.25;
rhow = 1027;
cd = 1.3;
a = 0.7859; b= 0.03477; c= 0.00412;
eps = 0.62197;
tau0 = 0.15;
K = 500;
k = Ly*K/H;
sbc = 5.670367 *10^(-8);
Ah = 10^(-3); % Take to be a constant but varies (Can be obtained from Table 2 of Smith 1988)
Pa = 10^5; 
alpha = 5.82*10^-5; 
beta = 8*10^-4; 


% define the grid
dz = 5;
dy =50*1000;

z = -[hm:dz:2500];
y = [0:dy:2000*1000];

% forcing conditions
tau = tau0*sin(pi*y/Ly);
Tair = (12-1)*y/Ly +1;
Simp = (35-34)*y/Ly + 34;
uw = sqrt(tau*rhow/rhoa/cd);
% end

%% Start algorithm to finally end up at the 





