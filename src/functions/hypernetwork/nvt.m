function [network, kappa_net] = nvt(network, FD, jamScenario, trscenario)
%NVT Summary of this function goes here
%   Detailed explanation goes here

% demandFactor = trscenario;
% x_nvt = [0.7456, 0.7707, 0.6329, 0.6560, 0.6594, 0.8201, 0.8813, 0.9470, 0.7871, 1];
% demandFactor = x_nvt(trscenario+1);
% demandFactor = 0.8;
demandFactor = 0.7707;
% demandFactor = 0.8;

% Define the numerical grid for a hyperlink
numGrid.dt = 0.1; % [s] Time-step. Spatial step is defined by dx=v_f*dt
numGrid.T = 7200; % [s] Simulation horizon:  Needs to be an integer of v_f
numGrid.numT = length(0:numGrid.dt:numGrid.T);
numGrid.dx = FD.u*numGrid.dt;
numGrid.precision = 10^-6;

links = network.links;

% Get corridor lengths
[group, id] = findgroups(links.corridor);
nCorrs = length(id);
func = @(length) sum(length);
corrLengths = splitapply(func, links.length, group);

hyperlink = struct();

for i = 1:nCorrs
    hyperlink(i).corridorLength = corrLengths(i);
    hyperlink(i).numT = numGrid.numT;
    hyperlink(i).numX = length(0:numGrid.dx:hyperlink(i).corridorLength);
    hyperlink(i).N = NaN(hyperlink(i).numX, numGrid.numT);
    hyperlink(i).BN = ones(hyperlink(i).numX, numGrid.numT)*FD.kc*FD.u;
    hyperlink(i).SF = zeros(hyperlink(i).numX, numGrid.numT);
    hyperlink(i).q = NaN(hyperlink(i).numX, numGrid.numT);
    hyperlink(i).k = NaN(hyperlink(i).numX, numGrid.numT);
    hyperlink(i).FD = FD;
end

% define nodes by connecting corridors, turning ratios, und pos on both
% corrs, offset, cycle length and red time
nNodes = max(network.connections.intersection);
idxNode = 1;
cycleTime = 90;
redTime = 45;

connections = network.connections;

for i = 1:22
    tempSubTable = connections(connections.intersection==i,:);
    nodes(idxNode).intersection = i;
    nodes(idxNode).corrs = unique(tempSubTable.fromCorr);
    actualTurnRatio = tempSubTable.Share(tempSubTable.fromCorr ~= tempSubTable.toCorr);
    if isempty(actualTurnRatio)
        nodes(idxNode).trs = [0;0;0;0];
    else
        nodes(idxNode).trs = actualTurnRatio;
    end
    
    fromLinks = unique(tempSubTable.fromLink);
    for j = 1:length(fromLinks)
        nodes(idxNode).positions(j) = links.cumLength(links.id==fromLinks(j));
    end
    nodes(idxNode).offsets = [0 45 0 45];
    nodes(idxNode).cycleTime = 90;
    nodes(idxNode).redTime = 45;
    nodes(idxNode).simTime = numGrid.T;
    nodes(idxNode).signal1 = VT_signal(cycleTime, redTime, nodes(idxNode).offsets(1), nodes(idxNode).positions(1), numGrid.T);
    nodes(idxNode).signal2 = VT_signal(cycleTime, redTime, nodes(idxNode).offsets(2), nodes(idxNode).positions(2), numGrid.T);
    nodes(idxNode).signal3 = VT_signal(cycleTime, redTime, nodes(idxNode).offsets(3), nodes(idxNode).positions(3), numGrid.T);
    nodes(idxNode).signal4 = VT_signal(cycleTime, redTime, nodes(idxNode).offsets(4), nodes(idxNode).positions(4), numGrid.T);
    idxNode = idxNode+1;
end

for i=1:nNodes
    hyperlink(nodes(i).corrs(1)) = putSignalsToBNs(numGrid, hyperlink(nodes(i).corrs(1)), nodes(i).signal1);
    hyperlink(nodes(i).corrs(2)) = putSignalsToBNs(numGrid, hyperlink(nodes(i).corrs(2)), nodes(i).signal2);
    hyperlink(nodes(i).corrs(3)) = putSignalsToBNs(numGrid, hyperlink(nodes(i).corrs(3)), nodes(i).signal3);
    hyperlink(nodes(i).corrs(4)) = putSignalsToBNs(numGrid, hyperlink(nodes(i).corrs(4)), nodes(i).signal4);
end

%% Definition of demand
% Definition of the initial density.

if jamScenario == 1
    demandFactor = 1;
end

for iCorr=1:nCorrs
    hyperlink(iCorr) = setInitialN(numGrid, hyperlink(iCorr), 0, [0 hyperlink(iCorr).corridorLength]);
    hyperlink(iCorr) = setUpstreamN(numGrid, hyperlink(iCorr), FD.kc*FD.u*0.5*demandFactor, [0 numGrid.T]);
end

%% Set the last BN to zero if there is a jam scenario
if jamScenario == 1
    for iCorr=1:nCorrs
        hyperlink(iCorr).BN(end-1,numGrid.T*0.75/numGrid.dt:end) = 0;
    end
end

%% Solve the given KWT problem
% Calculate N-Values for a hyperlink
[hyperlink] = solveVTnet(hyperlink, nodes, numGrid, FD);

% Calculate density and flows based on hyperlink
for iCorr=1:nCorrs
    hyperlink(iCorr) = calcDensFlow(numGrid, hyperlink(iCorr));
end

%% Extract data
[network, kappa_net] = extractData(FD, numGrid, hyperlink, network, jamScenario, nCorrs);

end

%% ===== Functions =====
function [network, kappa_net] = extractData(FD, numGrid, hyperlink, network, jamScenario, nCorrs)
kappa_net = NaN;
if jamScenario == 0
    for iCorr=1:nCorrs
        % Find location of BNs
        tmp = min(hyperlink(iCorr).BN, [], 2);
        tmp(end-1) = 0; % make sure the exit is also counted
        
        BN_temp = hyperlink(iCorr).BN(tmp==0,:);
        N_tmp = hyperlink(iCorr).N(tmp==0,:);
        N_checker = N_tmp(:,1:end-1) == N_tmp(:,2:end);
        
        sb_checker = N_checker & BN_temp(:,2:end) ~= 0; % BN(j) refers always to N(j-1)
        BN_temp2 = BN_temp;
        BN_temp2(BN_temp2==0) = NaN; 
      
        % Find flow reductions due to spillback only
        network.links.maxflow_i_sb(network.links.corridor == iCorr) = ...
            nanmean(BN_temp2(:,end-90/numGrid.dt:end-1) .* ...
            ~sb_checker(:,end-90/numGrid.dt:end-1),2) / FD.qmax;
        
        % record spillback times
        network.links.sb_i(network.links.corridor == iCorr) = ...
            sum(sb_checker(:,end-90/numGrid.dt:end-1),2);

        qTilde = (max(N_tmp(:,end-90/numGrid.dt:end-1),[],2)-min(N_tmp(:,end-90/numGrid.dt:end-1),[],2));
        network.links.maxflow_i(network.links.corridor == iCorr) = qTilde / FD.qmax ./ ...
            network.links.green(network.links.corridor == iCorr); % normalize it with FD and divide by green time
    end
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
else
    for iCorr=1:nCorrs
        k_corr(iCorr) = mean(hyperlink(iCorr).k(2:end-1,numGrid.T/numGrid.dt));
    end
    kappa_net = sum((k_corr .* [hyperlink.corridorLength]))/sum([hyperlink.corridorLength]);

end
end

function hyperlink = setInitialN(numGrid, hyperlink, densArray, posArray)
%SETINITIALN densArray is a nx1 array, posArray is a 2x1 array
% Convention: N(x=max, t=0) = 0
n_temp = -[0 cumsum(diff(posArray').*densArray')];
n_temp = n_temp + abs(min(n_temp));
for i =1:length(n_temp)-1
    N_temp = linspace(n_temp(i),n_temp(i+1), (diff(posArray(i,:)))/numGrid.dx+1);
    hyperlink.N(posArray(i,1)/numGrid.dx+1:posArray(i,2)/numGrid.dx,1)=N_temp(1:end-1)';
    hyperlink.N(end,1)=hyperlink.N(end-1,1) + (hyperlink.N(end-1,1) - hyperlink.N(end-2,1)); % Interpolation for very last value
end
end

function hyperlink = setUpstreamN(numGrid, hyperlink, flowArray, timeArray)
%SETUPSTREAMN flowArray is a nx1 array, timeArray is a 2x1 array
n_temp = [0 cumsum(diff(timeArray').*flowArray')];
for i =1:length(n_temp)-1
    N_temp = linspace(n_temp(i),n_temp(i+1), (diff(timeArray(i,:)))/numGrid.dt+1);
    hyperlink.N(1,timeArray(i,1)/numGrid.dt+1:timeArray(i,2)/numGrid.dt)=N_temp(1:end-1);
    hyperlink.N(1,end)=hyperlink.N(1,end-1) + (hyperlink.N(1,end-1) - hyperlink.N(1,end-2)); % Interpolation for very last value
end
end

function hyperlink = putSignalsToBNs(numGrid, hyperlink, signal)
%PUTSIGNALSTOBNS Save red lights to bottleneck matrix
signalsPosStep = signal.position./numGrid.dx+1;
tmpPhaseTimes = [0:numGrid.dt:numGrid.T]';
capacity = hyperlink.FD.kc*hyperlink.FD.u;

% Merge red and green
if min(signal.listGreenTime) > min(signal.listRedTime)
    all = [signal.listRedTime; signal.listGreenTime];
    tmp = all(:);
    tmpPhaseCosts = zeros(length(tmpPhaseTimes),1);
    for j=1:length(tmp)
        if mod(j,2)==1 % Is currently in a red phase
            tmpPhaseCosts(tmpPhaseTimes>tmp(j))=capacity;
        else  % Is currently in a green phase
            tmpPhaseCosts(tmpPhaseTimes>tmp(j))=0;
        end
    end
else
    all = [signal.listGreenTime; signal.listRedTime];
    tmp = all(:);
    tmpPhaseCosts = ones(length(tmpPhaseTimes),1)*capacity;
    for j=1:length(tmp)
        if mod(j,2)==1 % Is currently in a green phase
            tmpPhaseCosts(tmpPhaseTimes>tmp(j))=0;
        else  % Is currently in a red phase
            tmpPhaseCosts(tmpPhaseTimes>tmp(j))=capacity;
        end
    end
end

hyperlink.BN(signalsPosStep,:) = tmpPhaseCosts;
end

function hyperlink = calcDensFlow(numGrid, hyperlink)
%CALCDENSFLOW Calculates the flow and the density based on the N-values.
%   Fomrulas applied as in Leclercq & Paipuri (2019)

Ntemp = hyperlink.N;

% Flow
qtmp = [zeros(size(Ntemp,1),1) Ntemp(:,1:end-1)];
q = (Ntemp-qtmp)./numGrid.dt;

% Density
ktmp = [Ntemp(1,:); Ntemp(1:end-1,:)]; % assumptions density is constant in the very first x
k = -(Ntemp-ktmp)./numGrid.dx;

% feasibility constraints
q(q<0)=0;
q(q>hyperlink.FD.u*hyperlink.FD.kc)=hyperlink.FD.u*hyperlink.FD.kc;
k(k<0)=0;
k(k>hyperlink.FD.kappa)=hyperlink.FD.kappa;

% The code below is computationally more efficient than "round"
hyperlink.q = floor(q / numGrid.precision)*numGrid.precision;
hyperlink.k = round(k ,4); % this is needed for comparison
end



