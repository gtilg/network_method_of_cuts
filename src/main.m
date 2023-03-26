clear all;clc;

%% ----- (1) Define infrastructure input and initialize network -----
% Local FD parameters
FD.w = 5; % Wave speed [m/s]
FD.u = 10; % Free-flow speed [m/s]
FD.kappa = 0.15; % Jam density [veh/m]
FD.kc = FD.kappa * FD.w / (FD.w + FD.u); % critical density [veh/m]
FD.qmax = FD.kc * FD.u; % Capacity [veh/s]

% Load Network topology and control settings
tr_scenario = 1;
load('network_1.mat');

%% ----- (2) Create hypernetwork -----
% nVT approach
[network_nvt, ~] = nvt(network_sf_alt, FD, 0, tr_scenario);
[network_nvt, kappa_net] = nvt(network_nvt, FD, 1, tr_scenario);

% Approximate approach
approxLevel = 1;
[network_approx1] = create_hypernetwork_approx(network_sf_alt, FD, approxLevel, tr_scenario);
[kappa_net_app1, network_approx1] = findKappanet(network_approx1, FD);

approxLevel = 2;
[network_approx2] = create_hypernetwork_approx(network_sf_alt, FD, approxLevel, tr_scenario);
[kappa_net_app2, network_approx2] = findKappanet(network_approx2, FD);

%% ----- (3) Estimate the nMC MFD -----
% Estimate FF-branch, Estimate capacity branch, Estimate BW-branch
[tmp_mfd_nvt] = estimateFreeflowBranch(network_nvt, FD);
[tmp_mfd_nvt] = estimateCapacityBranch(tmp_mfd_nvt, network_nvt);
[mfd_net_nvt] = estimateCongestedBranch(tmp_mfd_nvt, kappa_net, FD);

[tmp_mfd_approx1] = estimateFreeflowBranch(network_approx1, FD);
[tmp_mfd_approx1] = estimateCapacityBranch(tmp_mfd_approx1, network_approx1);
[mfd_net_approx1] = estimateCongestedBranch(tmp_mfd_approx1, kappa_net_app1, FD);

[tmp_mfd_approx2] = estimateFreeflowBranch(network_approx2, FD);
[tmp_mfd_approx2] = estimateCapacityBranch(tmp_mfd_approx2, network_approx2);
[mfd_net_approx2] = estimateCongestedBranch(tmp_mfd_approx2, kappa_net_app2, FD);

%% ----- (4) Plotting -----
plot_method_comparison_all(mfd_net_nvt, mfd_net_approx1, mfd_net_approx2)
