function hyperlink = solveVTnet(hyperlink, nodes, numGrid, FD)
%SOLVEVT Solves a given VT problem
% The formula can be divided into three terms.
%   N(x,t) = min(N(x-Dx, t-Dt),N(x-Dx,t-psi*Dt)+kappaDx, N(x,t-Dt)+EtaDt)

kappa = FD.kappa;
dx = numGrid.dx;
dt = numGrid.dt;
nCorrs = length(hyperlink);

for i=1:nCorrs
    hpStruct(i).Ntemp = hyperlink(i).N;
    hpStruct(i).BNtemp = hyperlink(i).BN;
    hpStruct(i).SFtemp = hyperlink(i).SF;
    hpStruct(i).FeasMat = zeros(size(hyperlink(i).N));
end % Pause for-loop for the sake of clarity

%% Solve first two columns (t=dt)
for i=1:nCorrs
    hpStruct(i).Ntemp(2:end,2) = hpStruct(i).Ntemp(2:end,1);
end

%% Solve the rest
for j=3:numGrid.numT-1 % Go through time
    % Update values with results from previous timestep
    for i=1:nCorrs
        hpStruct(i).N_prev_t = hpStruct(i).Ntemp(:,j-1);
        hpStruct(i).N_prevprev_t = hpStruct(i).Ntemp(:,j-2);
    end
    
    % Do minimum operation
    hpStruct = calcNnextTimeStep(dx, dt, kappa, hpStruct, nodes, j);
    
    for i=1:nCorrs
        hpStruct(i).Ntemp(2:end-1,j)= hpStruct(i).N;
        % Find the feasibility matrix. It is used to exclude infeasbile
        % values for the flow and the density later on.
        hpStruct(i).FeasMat(:,j) = hpStruct(i).Ntemp(:,j)-hpStruct(i).Ntemp(:,j-1);
    end
end

% Store results
for i=1:nCorrs
    hyperlink(i).N = hpStruct(i).Ntemp;
end

%% Functions
    function hpStruct = calcNnextTimeStep(dx, dt, kappa, hpStruct, nodes, j)
        
        % Calculate general VT solution for all points in the numerical
        % grid
        % Links:
        [Ncand] = linkCalc(hpStruct, kappa, dx, dt, j);
        
        % Nodes: Adapt this solution at nodes with inflow/outflow
        [Ncand] = nodeCalc(Ncand, dx, kappa, hpStruct, nodes, j);
        
        % Apply the minimum operation
        for idx=1:length(hpStruct)
            tmp = [Ncand(idx).ff, Ncand(idx).bw, Ncand(idx).c,];
            hpStruct(idx).N = min(tmp(2:end-1,:),[],2);
        end
    end

    function [Ncand] = linkCalc(hpStruct, kappa, dx, dt, j)
        % Do VT calculation
        structSize = length(hpStruct);
        Ncand = struct('ff', cell(1,structSize), 'bw', cell(1, structSize), 'c', cell(1, structSize));
        for idx=1:length(hpStruct)
            % First term (free flow)
            Ncand(idx).ff = [NaN; hpStruct(idx).Ntemp(1:end-1,j-1)];
            % Second term (backward wave)
            Ncand(idx).bw = [hpStruct(idx).Ntemp(2:end,j-2) + kappa*dx; NaN];
            % Third term (stationary):
            Ncand(idx).c = hpStruct(idx).Ntemp(:,j-1) + hpStruct(idx).BNtemp(:,j).*dt;
        end
    end

    function [Ncand] = nodeCalc(Ncand, dx, kappa, hpStruct, nodes, j)
        % Loop over all nodes
        for jdx=1:length(nodes)
            % get relevant info from node properties
            cid=nodes(jdx).corrs;
            cpos=nodes(jdx).positions./dx+1; % +1 because first row is the boundary condition
            ctr=nodes(jdx).trs;
            
            % Free flow
            [Ncand] = nodeCalcff(Ncand, cid, cpos, ctr);
            
            % Backward wave
            if j > 3
                [Ncand] = nodeCalcbw(Ncand, hpStruct, cid, cpos, ctr, kappa, dx, j);
            end
        end
    end

    % Free flow case
    function [Ncand] = nodeCalcff(Ncand, cid, cpos, ctr)
        % This function considers inflows and outflows for downstream
        % sections.
        
        for idx=1:4
            N_main_end(idx) = Ncand(cid(idx)).ff(cpos(idx)+1) * (1 - ctr(idx));
            N_out(idx) = Ncand(cid(idx)).ff(cpos(idx)+1) * ctr(idx);
        end
        
        Ncand(cid(1)).ff(cpos(1)+1) = floor((N_main_end(1) + N_out(2)) / numGrid.precision) * numGrid.precision;
        Ncand(cid(2)).ff(cpos(2)+1) = floor((N_main_end(2) + N_out(3)) / numGrid.precision) * numGrid.precision;
        Ncand(cid(3)).ff(cpos(3)+1) = floor((N_main_end(3) + N_out(4)) / numGrid.precision) * numGrid.precision;
        Ncand(cid(4)).ff(cpos(4)+1) = floor((N_main_end(4) + N_out(1)) / numGrid.precision) * numGrid.precision;
    end

    % Backward wave case
    function [Ncand] = nodeCalcbw(Ncand, hpStruct, cid, cpos, ctr, kappa, dx, j)
        
        g = sum(hpStruct(cid(1)).BNtemp(cpos(1),j-2:j));
        g2 = hpStruct(cid(1)).BNtemp(cpos(1),j);
        
        for idx=1:4
            N_p_prime(idx) = hpStruct(cid(idx)).Ntemp(cpos(idx)+1,j-2); % downstream flow on corridor i, affected by congestion
            N_p_primeprime(idx) = hpStruct(cid(idx)).Ntemp(cpos(idx),j-3); % total free flow at ci
        end
        
        N_main = N_p_primeprime.*(1-ctr'); % what would have come in from ci in the free flow case
        
        % on corridor 1 flow would have come from corridor 2
        N_in(1) = hpStruct(cid(2)).Ntemp(cpos(2),j-3)*ctr(2);
        % on corridor 2 flow would have come from corridor 3
        N_in(2) = hpStruct(cid(3)).Ntemp(cpos(3),j-3)*ctr(3);
        % on corridor 3 flow would have come from corridor 4
        N_in(3) = hpStruct(cid(4)).Ntemp(cpos(4),j-3)*ctr(4);
        % on corridor 4 flow would have come from corridor 1
        N_in(4) = hpStruct(cid(1)).Ntemp(cpos(1),j-3)*ctr(1);
        
        delta = floor((N_p_prime - N_main - N_in) / numGrid.precision) * numGrid.precision;
        
        N_bw1_4 = N_p_primeprime(1) + 1/ctr(1)*(max([0, (dx*kappa+delta(4))*(0.5 - g2)/0.5, delta(4)+kappa*dx*(1-(1.5-g)/1.5*(1-ctr(4)))]));
        N_bw1_1 = N_p_primeprime(1) + 1/(1-ctr(1))*(max([0,(dx*kappa+delta(1))*(0.5 - g2)/0.5,delta(1)+kappa*dx*(1-(1.5-g)/1.5*ctr(4))]));
        
        N_bw2_1 = N_p_primeprime(2) + 1/ctr(2)*(max([0,(dx*kappa+delta(1))*(g2)/0.5, delta(1)+kappa*dx*(1-(g)/1.5*(1-ctr(1)))]));
        N_bw2_2 = N_p_primeprime(2) + 1/(1-ctr(2))*(max([0,(dx*kappa+delta(2))*(g2)/0.5,delta(2)+kappa*dx*(1-(g)/1.5*ctr(1))]));
        
        N_bw3_2 = N_p_primeprime(3) + 1/ctr(3)*(max([0, (dx*kappa+delta(2))*(0.5 - g2)/0.5, delta(2)+kappa*dx*(1-(1.5-g)/1.5*(1-ctr(2)))]));
        N_bw3_3 = N_p_primeprime(3) + 1/(1-ctr(3))*(max([0,(dx*kappa+delta(3))*(0.5 - g2)/0.5, delta(3)+kappa*dx*(1-(1.5-g)/1.5*ctr(2))]));
        
        N_bw4_3 = N_p_primeprime(4) + 1/ctr(4)*(max([0,(dx*kappa+delta(3))*(g2)/0.5, delta(3)+kappa*dx*(1-(g)/1.5*(1-ctr(3)))]));
        N_bw4_4 = N_p_primeprime(4) + 1/(1-ctr(4))*(max([0,(dx*kappa+delta(4))*(g2)/0.5,delta(4)+kappa*dx*(1-(g)/1.5*ctr(3))]));
        
        if ctr(1) == 1
            N_bw1_1 = N_p_primeprime(1) + kappa*dx*g/1.5;
        end
        
        if ctr(2) == 1
            N_bw2_2 = N_p_primeprime(2) + kappa*dx*(1.5-g)/1.5;
        end
        
        if ctr(3) == 1
            N_bw3_3 = N_p_primeprime(3) + kappa*dx*g/1.5;
        end
        
        if ctr(4) == 1
            N_bw4_4 = N_p_primeprime(4) + kappa*dx*(1.5-g)/1.5;
        end
        
        Ncand(cid(1)).bw(cpos(1)) = min(N_bw1_1, N_bw1_4);
        Ncand(cid(2)).bw(cpos(2)) = min(N_bw2_2, N_bw2_1);
        Ncand(cid(3)).bw(cpos(3)) = min(N_bw3_3, N_bw3_2);
        Ncand(cid(4)).bw(cpos(4)) = min(N_bw4_4, N_bw4_3);
    end
end
