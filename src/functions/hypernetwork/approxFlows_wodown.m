function [links] = approxFlows_wodown(x, links, nLinks, FD)
%APPROXSPILLBACKS Summary of this function goes here
%   Detailed explanation goes here

% Step 1: iterate to find demand
j = 1; % iterator
error = 1;
maxiter = 200;
for idxLink = 1:nLinks
    links.cumDemand(idxLink) = NaN; % preallocate
    cumDemand_new(idxLink) = NaN;
end
while (error >= 10^-5 || isnan(error)) && j <= maxiter
    for idxLink = 1:nLinks
        if links.origin(idxLink) == 1
            cumDemand_new(idxLink) = FD.qmax*x;
        else
            cumDemand_new(idxLink) = ...
                min(FD.qmax,...
                cumDemand_new(links.id==links.upstream_i(idxLink))*(1-links.alpha_ij(links.id==links.upstream_i(idxLink))) ...
                + cumDemand_new(links.id==links.upstream_j(idxLink))*links.alpha_ij(links.id==links.upstream_j(idxLink)));
        end
    end
    error = sum((abs([cumDemand_new]-[links.cumDemand]'))/abs([cumDemand_new]));
    j = j+1;
    for idxLink = 1:nLinks
        links.cumDemand(idxLink) = cumDemand_new(idxLink); % update the data
    end
end

% Step 2: Iterate to find spillbacks
j = 1; % iterator
error = 1;
maxiter = 100;

g = 45;
r = 45;

% Define standard values
boundaryCondition = struct;
spillbackData = struct;
for idxLink = 1:nLinks
    boundaryCondition(idxLink).q_downstream = FD.qmax;
    boundaryCondition(idxLink).q_up_i = FD.qmax;
    boundaryCondition(idxLink).q_up_j = FD.qmax;
    spillbackData(idxLink).sigma_idx = 0;
    spillbackData(idxLink).sigma_adjidx = 0;
    spillbackData(idxLink).flow = FD.qmax;
end

% Create temporary values for max flows
flow_temp = ones(1,nLinks)*FD.qmax;

while (error >= 10^-5 || isnan(error)) && j <= maxiter
    % Calculate spillbacks
    for idxLink = 1:nLinks
        % Calculate dN (as in the paper)
        [dN] = calcDN(idxLink, links, g, r, boundaryCondition, FD);
        % Calculate the spillback times based on dN
        [sigma_idx, sigma_adjidx] = calcSigmas(idxLink, links, dN);
        % Store the respective values
        spillbackData(idxLink).sigma_idx = sigma_idx;
        spillbackData(idxLink).sigma_adjidx = sigma_adjidx;
        spillbackData(idxLink).flow_i = max(0, FD.qmax / g * ( g - sigma_idx));
        spillbackData(idxLink).flow_j = max(0, FD.qmax / r * ( r - sigma_adjidx));
        
        % Roughly account for FIFO diverge
        if links.destination(idxLink) == 0
            spillbackData(idxLink).flow = min([spillbackData(idxLink).flow_i,spillbackData(links.id==links.j_link(idxLink)).flow_j]);
        else
            spillbackData(idxLink).flow = FD.qmax;
        end
    end
    
    % Check if change still occurs
    error = sum((abs(flow_temp-[spillbackData.flow]))/abs(flow_temp));
    flow_temp = [spillbackData.flow];
    
    % Update boundary condition
    for idxLink = 1:nLinks
        if links.destination(idxLink) == 1
            boundaryCondition(idxLink).q_downstream = FD.qmax;
            boundaryCondition(idxLink).q_up_i = FD.qmax;
            boundaryCondition(idxLink).q_up_j = FD.qmax;
        else
            % Do not update downstream propagation for this approximation
            % level
            % boundaryCondition(idxLink).q_downstream = spillbackData(links.id==links.downstream_i(idxLink)).flow;
            boundaryCondition(idxLink).q_up_i = spillbackData(idxLink).flow;
            boundaryCondition(idxLink).q_up_j = spillbackData(links.id==links.j_link(idxLink)).flow;
        end
    end
    
    j = j+1;
end

for idxLink = 1:nLinks
    links.q_sb(idxLink) = spillbackData(idxLink).flow;
end

q_sb = [spillbackData.flow];

% Step 3: iterate to find the minimum of demand and supply
j = 1; % iterator
error = 1;
maxiter = 100;
for idxLink = 1:nLinks
    links.maxflow_i(idxLink) = 1; % update the data
    q_app_tmp(idxLink) = 1;
end
while (error >= 10^-5 || isnan(error)) && j <= maxiter
    for idxLink = 1:nLinks
        if links.origin(idxLink) == 1
            q_app_tmp(idxLink) = min(q_sb(idxLink),links.cumDemand(idxLink));
        elseif links.destination(idxLink) == 1
            q_app_tmp(idxLink) = q_app_tmp(links.id==links.upstream_i(idxLink)) * ...
                (1-links.alpha_ij(links.id==links.upstream_i(idxLink))) + ...
                q_app_tmp(links.id==links.upstream_j(idxLink)) * links.alpha_ij(links.id==links.upstream_j(idxLink));
        else
            q_app_tmp(idxLink) = min(q_sb(idxLink), ...
                q_app_tmp(links.id==links.upstream_i(idxLink))*(1-links.alpha_ij(links.id==links.upstream_i(idxLink))) + ...
                q_app_tmp(links.id==links.upstream_j(idxLink)) * links.alpha_ij(links.id==links.upstream_j(idxLink)));
        end
    end
    error = sum((abs([q_app_tmp]-[links.maxflow_i]'))/abs([q_app_tmp]));
    j = j+1;
    for idxLink = 1:nLinks
        links.maxflow_i(idxLink) = q_app_tmp(idxLink); % update the data
        links.maxflow_i_sb(idxLink) = spillbackData(idxLink).flow;
        links.sb_i(idxLink) = spillbackData(idxLink).sigma_idx; 
    end
end

qapp = [links.maxflow_i];

output = mean(qapp)*(-1800);

end

function [dN] = calcDN(idxLink, links, g, r, boundaryCondition, FD)

in_red_phase = links.h(idxLink) == 0;

if links.destination(idxLink) == 1
    dN = 0;
else
    N_demand = links.cumDemand(idxLink)*(1-links.alpha_ij(idxLink))*g ...
        + links.cumDemand(links.id==links.j_link(idxLink))*links.alpha_ji(idxLink)*r;
    if ~isnan(links.downstream_i(idxLink)) && links.destination(links.id == links.downstream_i(idxLink)) == 1
        N_capacity = (g + r) * FD.qmax;
    else
        N_capacity = g * boundaryCondition(idxLink).q_downstream;
    end
    Delta_N = max([0, N_demand - N_capacity]);
    if in_red_phase
        dN = Delta_N - links.tau(idxLink) * links.cumDemand(links.id==links.j_link(idxLink)) * links.alpha_ji(idxLink);
    elseif ~in_red_phase
        dN = Delta_N - links.tau(idxLink) * links.cumDemand(idxLink) * (1-links.alpha_ij(idxLink));
    end
end

end

function [sigma_idx, sigma_adjidx] = calcSigmas(idxLink, links, dN)
in_red_phase = links.h(idxLink) == 0;

if links.destination(idxLink) == 1
    sigma_idx = 0;
    sigma_adjidx = 0; 
else
    if in_red_phase
        denominator = links.cumDemand(idxLink) * (1-links.alpha_ij(idxLink))*(1+sign(dN)) + ...
            links.cumDemand(links.id==links.j_link(idxLink)) * links.alpha_ji(idxLink)*(1-sign(dN));
        sigma_idx = max(0, 2*dN / denominator);
        sigma_adjidx = max(0,min(links.tau(idxLink), links.tau(idxLink) + 2*dN / denominator));
    elseif ~in_red_phase
        denominator = links.cumDemand(idxLink) * (1-links.alpha_ij(idxLink))*(1-sign(dN)) + ...
            links.cumDemand(links.id==links.j_link(idxLink)) * links.alpha_ji(idxLink)*(1+sign(dN));
        
        sigma_idx = max(0,min(links.tau(idxLink), links.tau(idxLink) + 2*dN / denominator));
        sigma_adjidx = max(0, 2*dN / denominator);
    end
end

end
