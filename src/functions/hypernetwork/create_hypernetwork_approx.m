function [network] = create_hypernetwork_approx(network, FD, approxLevel, scenario)
%CREATE_HYPERNETWORK_APPROX Summary of this function goes here
%   Detailed explanation goes here
%% Spillback parameter calculation
% calculate all parameters te, h, and tau for each link.
% assumes the same green time everywhere
nLinks = height(network.links);

for i=1:nLinks % This assumes links to be in order of corridor travel direction
    if network.links.destination(i) == 0
        network.links.te(i) = mod(network.links.red(i+1) + network.links.length(i+1)/FD.w,network.links.cycle(i)) - network.links.offset(i);
        if network.links.te(i) < 0
            network.links.te(i) = network.links.te(i) + network.links.cycle(i);
        end
    else
        network.links.te(i)=NaN;
    end
    network.links.tau(i) = mod(network.links.te(i),network.links.red(i));
    network.links.h(i) = max(network.links.te(i)+1-network.links.red(i),0)-max(network.links.te(i)-network.links.red(i),0); % zero if in red time
end

%% Uncomment for optimization
% lb = 0.5;
% ub = 1.0;
% x0 = 0.7;
% fun = @(x)approxFlows_opt(x, network.links, nLinks, FD);
% options = optimoptions('fmincon','Display','iter', 'Algorithm', 'active-set');
% x = fmincon(fun,x0,[],[],[],[],lb,ub,[], options);

% x_app1 = x;
% x_app2 = x;

% Results from optimization runs:
x_app1 = 1;
x_app2 = [0.8642, 0.8270, 0.8539, 0.7787, 0.8980, 0.8454, 0.9284, 0.8009, 0.8581, 0.8044];

if approxLevel == 1
    [links] = approxFlows_wodown(x_app1, network.links, nLinks, FD);
elseif approxLevel == 2
    x_app2_choice = x_app2(scenario+1);
    [links] = approxFlows(x_app2_choice, network.links, nLinks, FD);
end

network.links = links;

network.links.maxflow_i = network.links.maxflow_i/FD.qmax;
network.links.maxflow_i_sb = network.links.maxflow_i_sb/FD.qmax;
for j = 1:height(network.links)
    if network.links.destination(j) == 1
        network.links.maxflow_j(j) = network.links.maxflow_i(j);
        network.links.maxflow_j_sb(j) = network.links.maxflow_i_sb(j);
        network.links.sb_j(j) = network.links.sb_i(j);
    else
        network.links.maxflow_j(j) = network.links.maxflow_i(network.links.id == network.links.j_link(j));
        network.links.maxflow_j_sb(j) = network.links.maxflow_i_sb(network.links.id == network.links.j_link(j));
        network.links.sb_j(j) = network.links.sb_i(network.links.id == network.links.j_link(j));
    end
end

end

