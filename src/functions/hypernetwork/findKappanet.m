function [kappa_net_app, network] = findKappanet(network, FD)
%FINDKAPPANET Approximates kappa_net

tmp = network.links;

% 1. step: Find q=0 occurrence times at upstream and downstream end of
% links

% flip networks, as congestion travels backwards
network.links = flip(network.links);
tmp = flip(tmp);

% Assumes links to be in order of travel direction
sb_times = struct;
for i = 1:height(tmp)
    sb_times(i).id = network.links.id(i);
    sb_times(i).up = NaN;
    sb_times(i).down = NaN;
end

% with offsets, the first wave starts later (at the beginning of the red phase)
j = 1; % iterator
error = 1;
maxiter = 100;
while (error >= 10^-5 || isnan(error)) && j <= maxiter
    [sb_times_new] = findLinkSBtimes(sb_times, FD, tmp, network);
    error = sum((abs([sb_times_new.up]-[sb_times.up]))./abs([sb_times_new.up]));
    j = j+1;
    sb_times = sb_times_new; % update the data
end

% 2. Step: derive the according number of vehicles in the links. Assume a
% linear decrease with the time available to fill the link.
for i = 1:height(tmp)
    if tmp.origin(i) ~=1
        tmp.kappa_link(i) = FD.kappa * ...
            max([min([sb_times(i+1).supply_down, sb_times(i+1).demand_down])-(sb_times(i).down - tmp.length(i)/FD.u),0]) / ...
            (sb_times(i).supply_up-(sb_times(i).down - tmp.length(i)/FD.u));
    else
        tmp.kappa_link(i) = FD.kappa;
    end
end

% 3. step: Aggregate length weighted
kappa_net_app = sum([tmp.kappa_link].*[tmp.length])/sum([tmp.length]);
network.links.kappa_link = tmp.kappa_link;
network.links = flip(network.links);

end

function [sb_times] = findLinkSBtimes(sb_times, FD, tmp, network)
% (1) ====== supply-related constraints ======
for i = 1:height(tmp)
    if tmp.destination(i) == 1
        sb_times(i).supply_down = 0;
        % at the end demand is always decisive, thus:
        % "i+1" corresponds to furtehr downstream in this case
        sb_times(i).supply_up = tmp.length(i) * (FD.kappa-FD.kc) / (mean([tmp.maxflow_i(i+1),tmp.maxflow_j(i+1)]) * FD.qmax);
        sb_times(i).signal_down = 0;
    else
        % Find both downstream links and check for t_sb
        sb_down = min([sb_times([sb_times.id] == network.links.downstream_i(i)).up, ...
            sb_times([sb_times.id] == network.links.downstream_j(i)).up]);
        sb_times(i).supply_down = sb_down; % Store for later processing
        
        % Check if sb occurs during a red phase. If yes, the start of the red phase is the new start.
        [green, greenShare] = signalPhaseChecker(sb_down,network.links.offset(i),network.links.cycle(i),network.links.red(i));
        
        if green
            sb_times(i).signal_down = sb_down;
        else
            sb_times(i).signal_down = sb_down - mod(sb_down, network.links.cycle(i)) + network.links.offset(i);
        end
        
        % Check if a spillback existed during capacity state is occuring
        % now in the gridlock
        if tmp.origin(i) == 0 && tmp.sb_i(i+1) == 0 && tmp.sb_j(i+1) == 0
            % no spillbacks exist, demand is decisive
            sb_times(i).supply_up = sb_times(i).signal_down + tmp.length(i) * ...
                (FD.kappa - mean([tmp.maxflow_i(i+1),tmp.maxflow_j(i+1)]) * FD.qmax / FD.u)  / (mean([tmp.maxflow_i(i+1),tmp.maxflow_j(i+1)]) * FD.qmax);
        elseif tmp.origin(i) == 1
            sb_times(i).supply_up = sb_times(i).signal_down + tmp.length(i) * ...
                (FD.kappa - FD.qmax / FD.u)  / (FD.qmax);
        else
            % spillback is decisive
            sb_times(i).supply_up = sb_times(i).signal_down + tmp.length(i)/FD.w + greenShare * (tmp.red(i) - network.links.sb_i(i+1) - network.links.sb_j(i+1));
        end
    end
end

% (2) ====== demand-related constraints ======
for i = 1:height(tmp)
    if tmp.origin(i) == 1 || tmp.destination(i) == 1
        % Origin links get always fully congested, and destinations links
        % cannot propagate those effects further downstream, because they
        % are the most downstream
        sb_times(i).demand_down = NaN;
        sb_times(i).demand_up = NaN;
    else
        if i == 92
            1;
        end
        % linear interpoalte remaining vehicles into green phase
        time_rem_rel = mod(sb_times(i+1).supply_down + tmp.length(i)/FD.u, network.links.cycle(i));
        time_rem_rel_red = min(time_rem_rel,sb_times(i+1).supply_down - sb_times(i+1).signal_down);
        flow_rem_rel = tmp.maxflow_i(i+1) * (time_rem_rel - time_rem_rel_red)/tmp.green(i+1) + tmp.maxflow_j(i+1) * time_rem_rel_red/tmp.red(i+1);
        sb_times(i).demand_down = (time_rem_rel * flow_rem_rel)/ network.links.green(i)...
            + network.links.red(i) + floor((sb_times(i+1).supply_down + tmp.length(i)/FD.u)/network.links.cycle(i))*network.links.cycle(i);
        sb_times(i).demand_up = sb_times(i+1).demand_down;
    end
end

% (3) ====== Find constraints ======
for i = 1:height(tmp)
    sb_times(i).down = min([sb_times(i).supply_down, sb_times(i).demand_down, sb_times(i).signal_down]);
    sb_times(i).up = min([sb_times(i).supply_up, sb_times(i).demand_up]);
end
end

function [green, greenShare] = signalPhaseChecker(t,offset,cycle,red)
t_cycleSpecific = mod(t,cycle) - offset;
if t_cycleSpecific < 0
    t_cycleSpecific = t_cycleSpecific + cycle;
end
greenShare = max(0,(t_cycleSpecific - red)/45);
green = logical(max(t_cycleSpecific+1 - red,0)-max(t_cycleSpecific-red,0)); % green = 1 --> t in green phase
end


