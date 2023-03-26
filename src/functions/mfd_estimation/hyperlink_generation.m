function [ CalculationMethod, FD, HyperlinkProperties, Link ] = hyperlink_generation( CalculationMethod, FD, HyperlinkProperties, links, iCorr )
%HYPERLINK_GENERATION automatically completes the inputs defined by the
%user, and define the links of the hyperlink
[Link, HyperlinkProperties] = complete_Link_parameters( HyperlinkProperties, links, iCorr); % define the link properties
% ----- Link structure -----
% Link.ID : Link ID
% Link.Length : Link length [m]
% Link.Position : Global position of traffic signal [m] (auto)
% Link.NumLanes : Number of Lanes
% Link.SignalGT : Green Time [s]
% Link.SignalC  : Signal cycle [s]
% Link.SignalDelay : lag for the signal becoming green inside the cycle [s]
% Link.ListRedTime : List of times when the traffic signal turns red [s] (auto)
% Link.ListGreenTime : List of times when the traffic signal turns green [s] (auto)
% Link.PublicTransport : 1 if PT is active on link
% Link.ListPublicTransportStart : List of times when pt bottleneck
% starts
% Link.ListPublicTransportEnd : List of times when pt bottleneck
% ends

LastLink.Length = 250;
HyperlinkProperties.TotalLength = sum([Link.Length]) + LastLink.Length; % we add this information
HyperlinkProperties.G_Cycle = mean([Link(:).SignalC]); % we add this information
HyperlinkProperties.G_Ratio = mean([Link(:).SignalGT])/mean([Link(:).SignalC]); % we add this information
HyperlinkProperties.G_Offset = mean([Link(:).SignalDelay]); % we add this information
% ----- Hyperlink structure -----
% HyperlinkProperties.NumLink : number of link with signal at their end
% HyperlinkProperties.NumLane : number of lane on each link
% HyperlinkProperties.Length : average length of links [m]
% HyperlinkProperties.G_Cycle : common cycle for each signal [s]
% HyperlinkProperties.G_Ratio : average ratio green time / cycle for each signal
% HyperlinkProperties.G_Offset : offset for regular hyperlink [s]
% HyperlinkProperties.TypeChoice : 'Gambetta', 'OneCycleRand', 'Regular', or 'Specific'
% HyperlinkProperties.TotalLength : total length of the hyperlink [m]

CalculationMethod = adjust_time_window (CalculationMethod, FD, HyperlinkProperties);
% ----- CalculationMethod structure -----
% CalculationMethod.TimeWindow : duration of definition of signals [s]
% CalculationMethod.InflowType : form of the demand
% CalculationMethod.Error : value to avoid numerical errors

Link = automatic_settings( CalculationMethod, Link, LastLink, HyperlinkProperties );
end

function [Link, HyperlinkProperties] = complete_Link_parameters( HyperlinkProperties, links, iCorr )
%COMPLETE_LINK_PARAMETERS Puts link information in right structure format
%   Detailed explanation goes here

[links_in_iCorr] = find([links.corridor] == iCorr);
HyperlinkProperties.NumLink = length(links_in_iCorr); % number of links
Link = struct;
for i=1:HyperlinkProperties.NumLink
    index = links_in_iCorr(i);
    Link(i).ID = i; % Link Id
    Link(i).Length = links.length(index);
    Link(i).SignalGT = links.green(index);
    Link(i).SignalDelay = links.offset(index);
    Link(i).SignalC = links.cycle(index);
    Link(i).maxFlowIn = links.maxflow_j_sb(index); % max flow which comes in from other corridor, between 1 and 0
    Link(i).maxFlowOut = links.maxflow_i_sb(index); % max flow which would leave the current corridor, between 1 and 0
    Link(i).trIn = links.alpha_ji(index); % turning ratio which comes from other corridor, between 1 and 0
    Link(i).trOut = links.alpha_ij(index); % turning ratio which leaves the current corridor, between 1 and 0
    Link(i).spillbackStart = 9000;
    Link(i).spillbackEnd = 9000;
    
    if i==1
        Link(i).Position = Link(i).Length; % Global position of traffic signal [m]
    else
        Link(i).Position = Link(i-1).Position + Link(i).Length;
    end
    
    Link(i).NumLanes = HyperlinkProperties.NumLane; % Number of Lanes
    
    if HyperlinkProperties.TypeChoice.PublicTransport
        Link(i).PublicTransport = 1;
    else
        Link(i).PublicTransport = 0;
    end
end

end

function Link = automatic_settings( CalculationMethod, Link, LastLink, HyperlinkProperties )
%AUTOMATIC_SETTINGS automatically computes the other link variables
%(Link position, List of green and red times) based on previously defined
%variables

% ----- Completion of the Link parameters -----
% Create the lists for horiztonal nodes in the global VG
for i=1:1:HyperlinkProperties.NumLink
    
    % Signal times
    Link(i).ListGreenTime = Link(i).SignalDelay-(Link(i).SignalC*2) : Link(i).SignalC : CalculationMethod.TimeWindow; % List of times when the traffic signal turns green [s]
    Link(i).ListRedTime = Link(i).SignalDelay-(Link(i).SignalC*2)+Link(i).SignalGT : Link(i).SignalC : CalculationMethod.TimeWindow; % List of times when the traffic signal turns red [s]
    
    % Spillback bottlenecks
    Link(i).ListSpillbackStart = Link(i).spillbackStart-(2*Link(i).SignalC) : Link(i).SignalC : CalculationMethod.TimeWindow;  % List of times when the spillback begins to block [s]
    Link(i).ListSpillbackEnd = Link(i).spillbackEnd-(2*Link(i).SignalC) : Link(i).SignalC : CalculationMethod.TimeWindow; % List of times when the spillback ends to block [s]
    
    % PT bottlenecks
    if i == 1
        Link(i).ListPublicTransportStart = HyperlinkProperties.PT.start : HyperlinkProperties.PT.headway : CalculationMethod.TimeWindow;
        Link(i).ListPublicTransportEnd = HyperlinkProperties.PT.start + Link(i).Length/HyperlinkProperties.PT.speed : HyperlinkProperties.PT.headway : CalculationMethod.TimeWindow;
    else
        Link(i).ListPublicTransportStart = Link(i-1).ListPublicTransportEnd - Link(i).Length/HyperlinkProperties.PT.speed : HyperlinkProperties.PT.headway : CalculationMethod.TimeWindow;
        Link(i).ListPublicTransportEnd = Link(i-1).ListPublicTransportEnd + Link(i).Length/HyperlinkProperties.PT.speed : HyperlinkProperties.PT.headway : CalculationMethod.TimeWindow;
    end

    [hEdgeList] = merge_lists( Link(i), HyperlinkProperties );
    Link(i).ListHorizontalEdges = hEdgeList;
end

% --- last link (from upstream to downstream)
i = HyperlinkProperties.NumLink + 1;
Link(i).ID = i; % Link Id
Link(i).Length = LastLink.Length; % Link length [m]
Link(i).NumLanes = HyperlinkProperties.NumLane; % Number of Lanes
Link(i).SignalC = CalculationMethod.TimeWindow; % Signal cycle [s]
Link(i).SignalGT = CalculationMethod.TimeWindow; % Green Time [s]
Link(i).SignalDelay = 0; % lag for the signal becoming green inside the cycle [s]
Link(i).Position = Link(i-1).Position + Link(i).Length; % Position of the route end [m]
Link(i).maxFlowIn = Link(i-1).maxFlowIn;
Link(i).maxFlowOut = Link(i-1).maxFlowOut;
Link(i).trIn = 0;
Link(i).trOut = 0;
Link(i).spillbackStart = CalculationMethod.TimeWindow;
Link(i).spillbackEnd = CalculationMethod.TimeWindow;
Link(i).ListGreenTime = 0; % This theoretical last signal is always green [s]
Link(i).ListRedTime = CalculationMethod.TimeWindow; % This theoretical last signal is always green [s]
Link(i).ListPublicTransportStart = CalculationMethod.TimeWindow; % no bus on last theoretical link
Link(i).ListPublicTransportEnd = CalculationMethod.TimeWindow;
Link(i).ListSpillbackStart = CalculationMethod.TimeWindow; % no spillbacks on last theoretical link
Link(i).ListSpillbackEnd = CalculationMethod.TimeWindow;
Link(i).ListHorizontalEdges = CalculationMethod.TimeWindow;
Link(i).PublicTransport = 0;
end

function CalculationMethod = adjust_time_window (CalculationMethod, FD, HyperlinkProperties)
%ADJUST_TIME_WINDOW allows to be sure that the MFD will be correctly
%estimated

NumLinks = HyperlinkProperties.NumLink + 1; % link number
Cycle = HyperlinkProperties.G_Cycle; % hyperlink cycle
TotalLength = HyperlinkProperties.TotalLength; % hyperlink length

RequiredTimeWindow = TotalLength / FD.w + (NumLinks-1) * Cycle; % in the worst case, we are sure that platoons can be computed
%RequiredTimeWindow = CalculationMethod.TimeWindow;
CalculationMethod.TimeWindow = max(CalculationMethod.TimeWindow, RequiredTimeWindow); % we modify, if needed, the time window
end

function [hEdgeTimesList] = merge_lists( Link, HyperlinkProperties )
gCode = 1; rCode = 2; sCode = 3; ptCode = 4;
gCodeList = ones(1,length(Link.ListGreenTime))*gCode;
rCodeList = ones(1,length(Link.ListRedTime))*rCode;
sCodeListStart = ones(1,length(Link.ListSpillbackStart))*sCode;
sCodeListEnd = ones(1,length(Link.ListSpillbackEnd))*sCode;
ptCodeListStart = ones(1,length(Link.ListPublicTransportStart))*ptCode;
ptCodeListEnd = ones(1,length(Link.ListPublicTransportEnd))*ptCode;

if HyperlinkProperties.TypeChoice.PublicTransport == 1
    allTimesList = [Link.ListGreenTime Link.ListRedTime Link.ListSpillbackStart Link.ListSpillbackEnd Link.ListPublicTransportStart Link.ListPublicTransportEnd; ...
        gCodeList rCodeList sCodeListStart sCodeListEnd ptCodeListStart ptCodeListEnd];
elseif HyperlinkProperties.TypeChoice.PublicTransport == 0
    allTimesList = [Link.ListGreenTime Link.ListRedTime Link.ListSpillbackStart Link.ListSpillbackEnd; ...
        gCodeList rCodeList sCodeListStart sCodeListEnd];
end

[~, order] = sort(allTimesList(1,:));
hEdgeTimesList = allTimesList(:,order);
keep=ones(1,length(hEdgeTimesList));
for i=1:length(hEdgeTimesList)-1
    if hEdgeTimesList(1,i) == hEdgeTimesList(1,i+1) && hEdgeTimesList(2,i) ~= hEdgeTimesList(2,i+1) % e.g. spillback ends at red time.
        keep(i+1)=0;
    elseif hEdgeTimesList(1,i) == hEdgeTimesList(1,i+1) && hEdgeTimesList(2,i) == hEdgeTimesList(2,i+1) % e.g. spillback start == spillback end
        keep(i)=0;
        keep(i+1)=0;
    end
end
hEdgeTimesList=hEdgeTimesList(:,logical(keep));
end