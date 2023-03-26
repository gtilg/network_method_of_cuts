function [mfd] = estimateFreeflowBranch(network, FD)
%ESTIMATEFREEFLOWBRANCH Estimates the free flow branch of the MFD
% This code is partly based on the code corresponding to the paper 
% Leclercq, L., & Geroliminis, N. (2013). Estimating MFDs in simple 
% networks with route choice. 
% Transportation Research Part B: Methodological, 57, 468-484. Please cite
% this paper if you use any of the codes below.

hyperlinks = network.links;

%% ===== (1) General parameters =====
% ----- (0) Estimation parameters for the MFD and TSD Computation -----
CalculationMethod.TimeWindow = 3000; % please set a value high enough. Otherwise, it is the main cause of bug
CalculationMethod.BoundaryPaths = 1; % 1 = we consider the boundary paths on the first and last links, 0 otherwise (classical method)
CalculationMethod.EndRed = 0; % 1 = we only consider paths between two end of red phases, 0 = we consider all possible paths
% It is incompatible with the option "BoundaryPaths = 1"
CalculationMethod.Error = 10^-5; % admissible numerical error. Do not change this value...

% ----- (2) Global hyperlink parameters -----
% Define the scenario parameters for a regular hyperlink
HyperlinkProperties.NumLane = 1; % number of lanes

% Define parameters for PT
HyperlinkProperties.TypeChoice.PublicTransport = 0; %if PT is existing
HyperlinkProperties.PT.headway = 5*60; % Average headway in seconds
HyperlinkProperties.PT.speed = 12.5*0.8; % Average speed of buses in m/s
HyperlinkProperties.PT.start = 10; % Arrival of the first bus

%% ===== Derive the cuts for each corridor =====
nCorrs = max([hyperlinks.corridor]);
for iCorr = 1:nCorrs
    % ----- Hyperlink generation -----
    [CalculationMethod, FD, HyperlinkProperties, Link] = hyperlink_generation(CalculationMethod, FD, HyperlinkProperties, hyperlinks, iCorr);
    
    %% ===== (2) Models computation and graphics =====
    % ----- (1) Generation of the variational network and MFD computation -----
    [~, ~, ~, ~, MacroFundDiag, ~] = MFD(CalculationMethod, FD, Link, HyperlinkProperties);
    [q,k] = deriveFreeflowBranchFromCuts(MacroFundDiag, FD);
    qall(iCorr,:)=q;
end

%% ===== Aggregate the corridor cuts to network ones =====
% Lets aggregate over the whole free flow branch, but with only spillback
% impacted costs
mfd.k = k;
mfd.q = mean(qall,1);

function [q,k] = deriveFreeflowBranchFromCuts(MacroFundDiag, FD)
    temp_k_cand = MacroFundDiag.k*1000;
    temp_q_cand = MacroFundDiag.q*3600;
    
    k = 1:1:FD.kappa*1000;
    q = NaN(size(k));
    
    % Create base matrix
    for i = 1:k(end)
        
        k_high = min(temp_k_cand(round(temp_k_cand,4)>=i));
        k_low = max(temp_k_cand(temp_k_cand<i));
        try
            q_high = temp_q_cand(temp_k_cand==k_high);
        catch
            q(i)=0;
        end
        q_low = temp_q_cand(temp_k_cand==k_low);
        if k_high - k_low >0
            slope= (q_high-q_low)/(k_high-k_low);
            q(i) = slope*(i-k_low)+q_low;
        elseif k_high == k_low
            q(i) = q_low;
        end
        
        if isnan(q(i))
            q(i) = 0;
        end
    end
end

end

