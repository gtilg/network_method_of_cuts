function [network] = create_hypernetwork_approx(network, FD, method)
%CREATE_HYPERNETWORK_APPROX Here, the network representation for the
% approximate methods is created.

%% Spillback parameter calculation
% Calculate all parameters te, h, and tau for each link.
% Assumes the same green time everywhere.
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

%% Find solution
% Results from optimization runs:
x_app1 = 1;
x_app2 = 0.8270;

% LS approach
if method == "LS"
    [links] = approx(x_app1, network.links, nLinks, FD, method);

% FS approach
elseif method == "FS"
    [links] = approx(x_app2, network.links, nLinks, FD, method);
end

% Calculate and store results
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

