function [ Link, Arc, Node, Cut, MacroFundDiag, ListLowerCuts ] = MFD( CalculationMethod, FD, Link, HyperlinkProperties)

% ===== Generation of the variational network =====
% --- Generation of the structure Arcs ---
Arc = arc_initialisation( Link, FD, CalculationMethod ); % at this time, Arc is the set of phases of all signals

% --- Generation of the structure Nodes and completion of the structure Arcs ---
[Node, Arc] = node_initialisation (Arc); % at this time, Node is the end of phases of all signals and are linked to associated arc_id.

% --- Completion of the graph by arcs at speeds u and w ---
[Node, Arc] = variational_network( Link, Arc, Node, FD, CalculationMethod );

for i=1:length(Arc)
    Arc(i).duration=[Arc(i).Paval(1) - Arc(i).Pamont(1)]; % Associated link ID for horizontal arcs
end

% [p, G] = plot_graph(Arc, Node);
[G] = plot_graph(Arc, Node);

% ===== MFD Computing =====
% --- Computing the shortest paths in the graph ---
Cut = shortest_variational_paths( Link, Arc, Node, CalculationMethod, G );

% --- Computing the tight cuts and so the MFD ---
[MacroFundDiag, ListLowerCuts] = mfd_by_lower_cuts( Cut );
% disp(['Computation of the MFD : ',num2str(toc),' sec']);tic;
end



function [ Arc ] = arc_initialisation( Link, FD, CalculationMethod )
% Link structure : defined the link main characteristics
    % Link.ID : Link ID
    % Link.Length : Link length [m]
    % Link.Position : Global position of traffic signal [m] (auto)
    % Link.NumLanes : Number of Lanes
    % Link.SignalGT : Green Time [s]
    % Link.SignalC  : Signal cycle [s]
    % Link.SignalDelay : lag for the signal becoming green inside the cycle [s]
    % Link.ListRedTime : List of times when the traffic signal turns red [s] (auto)
    % Link.ListGreenTime : List of times when the traffic signal turns green [s] (auto)
    % Link.Capacity : Downstream capacity due to traffic signal [veh/s] (auto)
% Arc structure : Defined the arc main characteristics
    % Arc.Pamont : Coordinates of the upstream point of the arc
    % Arc.Paval : Coordinates of the downstream point of the arc
    % Arc.RankX : Associated link ID for horizontal arcs
    % Arc.CostRate : CostRate associated to the arc (The total cost is equal to the costrate multiplied by the time distance of the arc
    % Arc.UsptreamNode : Node ID of the upstream point
    % Arc.DownstreamNode : Node ID of the downstream point

IdArc=1; % initialisation of arc index
Arc=struct; % initialisation of the Arc structure: delete the "orange-error" on the right side
NumLinks = length(Link) - 1; % Link number without counting the last virtual signal

for i=1:1:NumLinks
    % ---- for each link (from upstream to downstream)
    for j=1:1:length(Link(i).ListHorizontalEdges)
        
        if j < length(Link(i).ListHorizontalEdges)
            Arc(IdArc).Pamont=[Link(i).ListHorizontalEdges(1,j) Link(i).Position]; % Coordinates of the upstream point of the arc
            Arc(IdArc).Paval=[Link(i).ListHorizontalEdges(1,j+1) Link(i).Position]; % Coordinates of the downstream point of the arc
        else
            Arc(IdArc).Pamont=[Link(i).ListHorizontalEdges(1,j) Link(i).Position]; % Coordinates of the upstream point of the arc
            Arc(IdArc).Paval=[CalculationMethod.TimeWindow Link(i).Position]; % Coordinates of the downstream point of the arc
        end

        Arc(IdArc).RankX=i; % Associated link ID for horizontal arcs
        Arc(IdArc).duration=[Arc(IdArc).Paval(1) - Arc(IdArc).Pamont(1)]; % Associated link ID for horizontal arcs
        
        % Derive costs based on code
        edgeCode = Link(i).ListHorizontalEdges(2,j);
        if edgeCode == 1 % green edge follows
            Arc(IdArc).edgeCode = edgeCode;
            Arc(IdArc).CostRate=Link(i).NumLanes*FD.qmax*Link(i).maxFlowOut; % CostRate associated to the arc
        elseif edgeCode == 2 % red edge follows
            Arc(IdArc).edgeCode = edgeCode;
            Arc(IdArc).CostRate=0; % CostRate for a red phase
        elseif edgeCode == 3 % spillback edge follows
            Arc(IdArc).edgeCode = edgeCode;
            Arc(IdArc).CostRate=0; % CostRate for a red phase
        elseif edgeCode == 4 % PT edge follows
            Arc(IdArc).edgeCode = edgeCode;
            Arc(IdArc).CostRate=Link(i).NumLanes*0.5*FD.qmax; % CostRate for a red phase
        end
        
        IdArc=IdArc+1; % completed computation of this arc
    end
end % completed arc creation
end

function [Node, Arc] = node_initialisation (Arc)
% Node identification and connection to arcs:
    % Node.T : T-Coordinate of the nodes (t,x) plane
    % Node.X : X-Coordinate of the nodes (t,x) plane
    % Node.UpstreamArc : List of Id of the upstream arcs
    % Node.NumUpstreamArc : Number of upstream arcs
    % Node.DownstreamArc : List of Id of the downstream arcs
    % Node.NumDownstreamArc : Number of donwstream arcs
    % Node.Type:
        % 0, the node is connected with a downstream red phase;
        % 1, the node is connected with a downstream green phase;
        % 2 the node is a intermediate node created by intersecting oblique arc and will be the origin of a wave at speed u
        % 3 the node is a intermediate node created by intersecting oblique arc and will be the origin of a wave at speed w
        % TODO: THAT NEEDS TO BE DOUBLE CHECKED
% Arc structure : Defined the arc main characteristics
    % Arc.Pamont : Coordinates of the upstream point of the arc
    % Arc.Paval : Coordinates of the downstream point of the arc
    % Arc.RankX : Associated link ID for horizontal arcs
    % Arc.CostRate : CostRate associated to the arc (The total cost is equal to the costrate multiplied by the time distance of the arc
    % Arc.UsptreamNode : Node ID of the upstream point
    % Arc.DownstreamNode : Node ID of the downstream point
   
% ----- Initialisation : P is a struct array containing all nodes in the network -----
P(1,:)=Arc(1).Pamont; % first node: the upstream point of the first arc
P(2,:)=Arc(1).Paval; % second node: the downstream point of the first arc
IdNode=3; % next node to define
NumArcs = length(Arc); % arc number at initialisation step
PDownArc = zeros(NumArcs,1); % Index of the arc down a node
PUpArc = zeros(NumArcs,1); % Index of the arc up a node

for i=1:1:NumArcs
    % ---- for each arc defined before
    
    % --- Early point of the considered Arc    
    IdMatch = find((P(:,1) == Arc(i).Pamont(1) & P(:,2) == Arc(i).Pamont(2))==1);
    % we search which node is the early point of the considered arc
    % (both coordinates) in the set of already defined nodes
    if isempty(IdMatch)
        % if the early point is not in the set of defined nodes, we add it
        P(IdNode,:) = Arc(i).Pamont; % addition of the early point
        IdMatch = IdNode; % we consider this new node
        IdNode = IdNode+1; % next node to define
    end
    PDownArc(IdMatch)=i; % the node indexed by "IdMatch" is the early point of the arc "i"
    
    % --- Late point of the considered Arc    
    IdMatch=find((P(:,1) == Arc(i).Paval(1) & P(:,2) == Arc(i).Paval(2))==1);
    % we search which node is the late point of the considered arc
    % (both coordinates) in the set of already defined nodes
    if isempty(IdMatch)
        % if the late point is not in the set of defined nodes, we add it
        P(IdNode,:)=Arc(i).Paval; % addition of the late point
        IdMatch=IdNode; % we consider this new node   
        IdNode=IdNode+1; % next node to define
    end    
    PUpArc(IdMatch)=i; % the node indexed by "IdMatch" is the late point of the arc "i"
end
% at the end of this loop, for each node, we found 1 or 0 arc up the node,
% and 1 or 0 arc down the node. Never more.

NumNodes=IdNode-1; % Node number
Node=struct; % initialisation of the Node structure: delete the "orange-error" on the right side



% ----- Creation of the node structure -----
for i=1:1:NumNodes
    % ---- for each node defined before
    Node(i).T = P(i,1); % T-Coordinate of the nodes (t,x) plane
    Node(i).X = P(i,2); % X-Coordinate of the nodes (t,x) plane
    
    % --- find an upstream arc
    if (i<=length(PUpArc))
        % if there is an arc up the considered node
        Node(i).UpstreamArc = PUpArc(i); % List of Id of the upstream arcs
    else
        % there is not any arc up the considered node
        Node(i).UpstreamArc = 0; % Convention
    end
    
    if Node(i).UpstreamArc==0
        Node(i).NumUpstreamArc = 0; % Number of upstream arcs
    else
        Node(i).NumUpstreamArc = 1; % Number of upstream arcs
        Node(i).RankX = Arc(Node(i).UpstreamArc).RankX; % Associated link ID
        Arc(PUpArc(i)).DownstreamNode = i; % Node ID of the downstream point
    end
    
    % --- find a downstream arc
    if (i<=length(PDownArc))
        % if there is an arc down the considered node
        Node(i).DownstreamArc = PDownArc(i); % List of Id of the downstream arcs
    else
        % there is not any arc up the considered node
        Node(i).DownstreamArc = 0; % List of Id of the downstream arcs
    end
    
    if Node(i).DownstreamArc==0
        Node(i).NumDownstreamArc = 0; % Number of donwstream arcs
    else
        Node(i).NumDownstreamArc = 1; % Number of donwstream arcs
        Node(i).RankX = Arc(Node(i).DownstreamArc).RankX;
        Arc(PDownArc(i)).UpstreamNode = i; % Node ID of the downstream point
    end
    
    % --- determination of the node type
    if Node(i).NumDownstreamArc > 0 && Arc(Node(i).DownstreamArc).CostRate>0
        Node(i).Type = 1; % 1, the node is connected with a downstream green phase;
    else
        Node(i).Type = 0; % 0, the node is connected with a downstream red phase;
        % at this stage, a node cannot be connected with both red and green phases
    end
end
end

function [ Node, Arc ] = variational_network( Link, Arc, Node, FD, CalculationMethod )
% Adding Oblique line into the network and connections

% Link structure : defined the link main characteristics
    % Link.ID : Link ID
    % Link.Length : Link length [m]
    % Link.Position : Global position of traffic signal [m] (auto)
    % Link.NumLanes : Number of Lanes
    % Link.SignalGT : Green Time [s]
    % Link.SignalC  : Signal cycle [s]
    % Link.SignalDelay : lag for the signal becoming green inside the cycle [s]
    % Link.ListRedTime : List of times when the traffic signal turns red [s] (auto)
    % Link.ListGreenTime : List of times when the traffic signal turns green [s] (auto)
    % Link.Capacity : Downstream capacity due to traffic signal [veh/s] (auto)

% Arc structure : Defined the arc main characteristics
    % Arc.Pamont : Coordinates of the upstream point of the arc
    % Arc.Paval : Coordianates of the downstream point of the arc
    % Arc.RankX : Associated link ID for horizontal arcs
    % Arc.CostRate : CostRate associated to the arc (The total cost is equal to the costrate multiplied by the time distance of the arc
    % Arc.UpstreamNode : Node ID of the upstream point
    % Arc.DownstreamNode : Node ID of the downstream point

% Node structure: Defined the node main characteristics and connections to arcs
    % Node.T : T-Coordinate of the nodes (t,x) plane
    % Node.X : X-oordinate of the nodes (t,x) plane
    % Node.UpstreamArc : List of Id of the upstream arcs
    % Node.NumUpstreamArc : Number of upstream arcs
    % Node.DownstreamArc : List of Id of the downstream arcs
    % Node.NumDownstreamArc : Number of donwstream arcs
    % Node.Type:
        % 0, the node is connected with a downstream red phase;
        % 1, the node is connected with a downstream green phase;
        % 2 the node is a intermediate node created by intersecting oblique arc and wille be the origin of a wave at speed u
        % 3 the node is a intermediate node created by intersecting oblique arc and wille be the origin of a wave at speed w

ListStartingNodeId= find(([Node(:).Type]==1) & [Node(:).T]>=0); % List of Starting points for the oblique line
% NB: Starting points correspond to the begining of new green phases, i.e. end of last red phases
% sorting the list of starting nodes from the earliest to the latest one
[SortedTime,SortedNodes] = sort([Node(ListStartingNodeId).T]);
ListStartingNodeId = ListStartingNodeId(SortedNodes);

IdList = 0; % next node to be considered as the start of new arcs (upstream at speed w and downstream at speed u)
TmpLengthList = length(ListStartingNodeId); % length of the list of nodes to be considerd
IdNewNode = length(Node) ; % next node to define
IdNewArc = length(Arc); % next arc to define
RankXMax = length(Link) - 1; % number of signals whithout counting the last virtual signal
TimeAllowance = 0.001; % [s] when we compare 

while (IdList < TmpLengthList)
    % while nodes ot be considered are remaining, we complete the variational network
    IdList = IdList + 1 ; % index of the investigated node
    CurrentNode = Node(ListStartingNodeId(IdList)); % Current Investigated node
    StartPoint = [CurrentNode.T CurrentNode.X]; % coordinates of the considered node
    
    % --- investigation downstream (at speed u) ---
    if (CurrentNode.RankX < RankXMax) && ((CurrentNode.Type==1) || (CurrentNode.Type==2))
        % for all signals exept for the last one
        NewNodeSpace = Link(CurrentNode.RankX+1).Position; % next upstream signal
        Temp_ObliqueSpeed = FD.u;
        Temp_ObliqueCostRate = 0;
        NewNodeTime = StartPoint(1)+(NewNodeSpace-StartPoint(2))/Temp_ObliqueSpeed; % time corresponding to the next point
        EndPoint = [NewNodeTime, NewNodeSpace]; % coordinates of the point downstream the arc
        if NewNodeTime<=CalculationMethod.TimeWindow-100 % QUICK AND DIRTY BUGFIX
            % then there is a downstream horizontal arc intersected by the wave u
            [Arc,Node,IdNewArc,IdNewNode,ListStartingNodeId,SortedTime,TmpLengthList] = complete_variational_network( StartPoint,EndPoint,Temp_ObliqueCostRate,Temp_ObliqueSpeed,Arc,IdNewArc,Node,CurrentNode,IdNewNode,NewNodeTime,NewNodeSpace,ListStartingNodeId,SortedTime,TmpLengthList,IdList,TimeAllowance );
        end
    end
    
    % --- investigation upstream (at speed w) ---
    if (CurrentNode.RankX > 1)  && ((CurrentNode.Type==1) || (CurrentNode.Type==3))
        % for all signals exept for the first one
        NewNodeSpace = Link(CurrentNode.RankX-1).Position; % next downstream signal
        Temp_ObliqueSpeed = -FD.w; % wave speed w
        Temp_ObliqueCostRate = Link(CurrentNode.RankX).NumLanes * FD.w * FD.kappa;
        NewNodeTime = StartPoint(1)+(NewNodeSpace-StartPoint(2))/Temp_ObliqueSpeed; % time corresponding to the next point
        EndPoint = [NewNodeTime, NewNodeSpace]; % coordinates of the point downstream the arc
        if NewNodeTime<=CalculationMethod.TimeWindow-100 % QUICK AND DIRTY BUGFIX: When slanted arc would go close to calculationwindow, it could be that there is actually no horizontal arc anymore
            % then there is a downstream horizontal arc intersected by the wave w
            [Arc,Node,IdNewArc,IdNewNode,ListStartingNodeId,SortedTime,TmpLengthList] = complete_variational_network( StartPoint,EndPoint,Temp_ObliqueCostRate,Temp_ObliqueSpeed,Arc,IdNewArc,Node,CurrentNode,IdNewNode,NewNodeTime,NewNodeSpace,ListStartingNodeId,SortedTime,TmpLengthList,IdList,TimeAllowance );
        end
    end
end

% completion of the network with extremal arcs
if CalculationMethod.BoundaryPaths
    [ Arc, Node ] = complete_network_boundaries( Link, Arc, Node, FD );
end
end

function [ Arc,Node,IdNewArc,IdNewNode,ListStartingNodeId,SortedTime,TmpLengthList ] = complete_variational_network( StartPoint,EndPoint,Temp_ObliqueCostRate,Temp_ObliqueSpeed,Arc,IdNewArc,Node,CurrentNode,IdNewNode,NewNodeTime,NewNodeSpace,ListStartingNodeId,SortedTime,TmpLengthList,IdList,TimeAllowance )
%COMPLETE_VARIATIONAL_NETWORK completes the variational network with a new
%arc and makes the necessary connections

% New oblique arc creation and node connection
IdNewArc = IdNewArc + 1; % next arc to define
Arc(IdNewArc).Pamont = StartPoint;
Arc(IdNewArc).UpstreamNode = ListStartingNodeId(IdList);
Arc(IdNewArc).Paval = EndPoint;
Arc(IdNewArc).RankX = 0; % Convention: an oblique arc has a rankx equal to 0
Arc(IdNewArc).CostRate = Temp_ObliqueCostRate;
Node(ListStartingNodeId(IdList)).DownstreamArc(end+1) = IdNewArc; % this arc is downstream the considered node
Node(ListStartingNodeId(IdList)).NumDownstreamArc = Node(ListStartingNodeId(IdList)).NumDownstreamArc + 1; % there is a new arc downstream the considered node

Tmp_Pamont = [Arc(:).Pamont];Tmp_Pamont = Tmp_Pamont(1:2:end);
Tmp_Paval = [Arc(:).Paval];Tmp_Paval = Tmp_Paval(1:2:end);

if Temp_ObliqueSpeed > 0
    % --- wave at speed u (>0) ---
    IntersectedArc = find(([Arc(:).RankX] == CurrentNode.RankX+1) & (Tmp_Pamont(1,:) <= NewNodeTime) & (Tmp_Paval(1,:) >= NewNodeTime),1,'first'); % Arc intersected by the wave u
else
    % --- wave at speed w (<0) ---
    IntersectedArc = find(([Arc(:).RankX] == CurrentNode.RankX-1) & (Tmp_Pamont(1,:) <= NewNodeTime) & (Tmp_Paval(1,:) >= NewNodeTime),1,'first'); % Arc intersected by the wave w
end

if (Arc(IntersectedArc).CostRate==0)
    ArcRedPhase = 1; % the intersecting arc corresponds to a red phase
else
    ArcRedPhase = 0; % the intersecting arc does NOT correspond to a red phase
end

% Check if the intersection point is already a Node
IdNode = find( (abs([Node(:).T]-NewNodeTime)<=TimeAllowance) & ([Node(:).X]==NewNodeSpace) );
if isempty(IdNode)
    % This is a new point
    IdNewNode = IdNewNode + 1; % new point to add
    Node(IdNewNode).T = NewNodeTime;
    Node(IdNewNode).X = NewNodeSpace;
    Arc(IdNewArc).DownstreamNode = IdNewNode;

    % We need to cut the horizontal intersecting arcs into 2
    ArcToCut = Arc(IntersectedArc);
    % --- first part of the arc to cut: we keep the same id
    Arc(IntersectedArc).Paval = EndPoint; % Coordianates of the downstream point of the arc
    Arc(IntersectedArc).DownstreamNode = IdNewNode; % Node ID of the downstream point
    % The other variables do not change
    % --- second part of the arc to cut: we add a new arc
    IdNewArc = IdNewArc + 1; % next arc to define
    Arc(IdNewArc).Pamont = EndPoint; % Coordinates of the upstream point of the arc
    Arc(IdNewArc).Paval = ArcToCut.Paval; % Coordinates of the downstream point of the arc
    Arc(IdNewArc).RankX = ArcToCut.RankX; % associated link ID for horizontal arcs
    Arc(IdNewArc).CostRate = ArcToCut.CostRate; % associated to the arc (The total cost is equal to the costrate multiplied by the time distance of the arc
    Arc(IdNewArc).duration=[Arc(IdNewArc).Paval(1) - Arc(IdNewArc).Pamont(1)]; % Associated link ID for horizontal arcs
    Arc(IdNewArc).UpstreamNode = IdNewNode; % Node ID of the upstream point
    Arc(IdNewArc).DownstreamNode = ArcToCut.DownstreamNode; % Node ID of the downstream point

    % completing the information for the new node
    if Temp_ObliqueSpeed > 0
        % --- wave at speed u (>0) ---
        Node(IdNewNode).RankX = CurrentNode.RankX + 1;
    else
        % --- wave at speed w (<0) ---
        Node(IdNewNode).RankX = CurrentNode.RankX - 1;
    end
    Node(IdNewNode).UpstreamArc = [IntersectedArc, IdNewArc-1]; % List of Id of the upstream arcs: the first part of the intersected arc, and the oblique arc
    Node(IdNewNode).NumUpstreamArc = 2; % Number of upstream arcs: the first part of the intersected arc, and the oblique arc
    Node(IdNewNode).DownstreamArc = IdNewArc; % List of Id of the upstream arcs: the second part of the intersected arc
    Node(IdNewNode).NumDownstreamArc = 1; % Number of upstream arcs: the second part of the intersected arc
    if Temp_ObliqueSpeed > 0
        % --- wave at speed u (>0) ---
        Node(IdNewNode).Type = 2; % Convention
    else
        % --- wave at speed w (<0) ---
        Node(IdNewNode).Type = 3; % Convention
    end

    if ~ArcRedPhase
        % if the new point is during a green phase, it has to
        % be added to point to consider
        CutScanningNode = find(SortedTime(IdList:end)<NewNodeTime,1,'last'); % the first element corresponds
        if ~isempty(CutScanningNode)
            CutScanningNode = CutScanningNode + IdList - 1; % in the classical case
        else
            CutScanningNode = IdList; % in the case of infinite speed u
        end
        SortedTime = [SortedTime(1:CutScanningNode), NewNodeTime ,SortedTime(CutScanningNode+1:end)];
        ListStartingNodeId = [ListStartingNodeId(1:CutScanningNode), IdNewNode ,ListStartingNodeId(CutScanningNode+1:end)];
        TmpLengthList = TmpLengthList + 1;
    end
else
    % This point already exists: node update
    Arc(IdNewArc).DownstreamNode = IdNode;
    Node(IdNode).UpstreamArc(end+1) = IdNewArc; % the considered arc is added
    Node(IdNode).NumUpstreamArc = Node(IdNode).NumUpstreamArc + 1; % the considered arc is added
    Node(IdNode).Type = 1;
end
end

function [ Arc, Node ] = complete_network_boundaries( Link, Arc, Node, FD )
%COMPLETE_NETWORK_BOUNDARIES completes the variational network by adding the final
%arcs at speed u and w. Moreover, it adds extremal Nodes

NumArcs = length(Arc);
NumNodes = length(Node);
RankXMax = length(Link) - 1; % number of signals whithout counting the last virtual signal

% ----- (1) Completion of ascending paths on the first link -----
Condition1 = ([Node.RankX] == 1); % nodes on the first signal
Condition2 = ([Node.Type] == 1); % ends of red phases
Liste = find(Condition1 & Condition2);
for i_node = 1 : length(Liste)
    % for each node satisfaying the above conditions
    % we have to add an upstream arc (speed u) and an upstream node
    NumNodes = NumNodes + 1; % new node
    NumArcs = NumArcs + 1; % new arc
    % --- (a) we complete the information about the current node ---
    Node(Liste(i_node)).UpstreamArc(end+1) = NumArcs;
    Node(Liste(i_node)).NumUpstreamArc = Node(Liste(i_node)).NumUpstreamArc + 1;
    CurrentNode = Node(Liste(i_node));
    % --- (b) we add the final node ---
    Node(NumNodes).T = CurrentNode.T + (0 - CurrentNode.X) / FD.u; % time corresponding to the next point
    Node(NumNodes).X = 0; % position of the new node : start of the arterial
    Node(NumNodes).UpstreamArc = [];
    Node(NumNodes).NumUpstreamArc = 0;
    Node(NumNodes).DownstreamArc = NumArcs; % id of the new arc
    Node(NumNodes).NumDownstreamArc = 1;
    Node(NumNodes).RankX = 0; % convention : start of the arterial
    Node(NumNodes).Type = 1 ; % origin of a path (necessary for begin a path)
    NewNode = Node(NumNodes);
    % --- (c) we add the final arc ---
    Arc(NumArcs).Pamont = [NewNode.T NewNode.X];
    Arc(NumArcs).Paval = [CurrentNode.T CurrentNode.X];
    Arc(NumArcs).RankX = 0; % Convention: an oblique arc has a rankx equal to 0
    Arc(NumArcs).CostRate = 0; % wave at speed u
    Arc(NumArcs).duration=[Arc(NumArcs).Paval(1) - Arc(NumArcs).Pamont(1)];
    Arc(NumArcs).UpstreamNode = NumNodes;
    Arc(NumArcs).DownstreamNode = Liste(i_node);
end

% ----- (2) Completion of ascending paths on the last link -----
Condition1 = ([Node.RankX] == RankXMax); % nodes on the first signal
Condition2a = ([Node.Type] == 1); % ends of red phases
Condition2b = ([Node.Type] == 2); % ends of arcs at the speed u
Liste = find(Condition1 & (Condition2a | Condition2b));
for i_node = 1 : length(Liste)
    % for each node satisfaying the above conditions
    % we have to add a downstream arc (speed u) and a downstream node
    NumNodes = NumNodes + 1; % new node
    NumArcs = NumArcs + 1; % new arc
    % --- (a) we complete the information about the current node ---
    Node(Liste(i_node)).DownstreamArc(end+1) = NumArcs;
    Node(Liste(i_node)).NumDownstreamArc = Node(Liste(i_node)).NumDownstreamArc + 1;
    CurrentNode = Node(Liste(i_node));
    % --- (b) we add the final node ---
    Node(NumNodes).T = CurrentNode.T + (Link(end).Position - CurrentNode.X) / FD.u; % time corresponding to the next point
    Node(NumNodes).X = Link(end).Position; % position of the new node : start of the arterial
    Node(NumNodes).UpstreamArc = NumArcs; % id of the new arc
    Node(NumNodes).NumUpstreamArc = 1;
    Node(NumNodes).DownstreamArc = [];
    Node(NumNodes).NumDownstreamArc = 0;
    Node(NumNodes).RankX = RankXMax + 1; % convention : end of the arterial
    Node(NumNodes).Type = 2 ; % end of a wave at speed u
    NewNode = Node(NumNodes);
    % --- (c) we add the final arc ---
    Arc(NumArcs).Pamont = [CurrentNode.T CurrentNode.X];
    Arc(NumArcs).Paval = [NewNode.T NewNode.X];
    Arc(NumArcs).RankX = 0; % Convention: an oblique arc has a rankx equal to 0
    Arc(NumArcs).CostRate = 0; % wave at speed u
    Arc(NumArcs).duration=[Arc(NumArcs).Paval(1) - Arc(NumArcs).Pamont(1)];
    Arc(NumArcs).UpstreamNode = Liste(i_node);
    Arc(NumArcs).DownstreamNode = NumNodes;
end

% ----- (3) Completion of descending paths on the last link -----
Condition1 = ([Node.RankX] == RankXMax); % nodes on the first signal
Condition2 = ([Node.Type] == 1); % ends of red phases
Liste = find(Condition1 & Condition2);
for i_node = 1 : length(Liste)
    % for each node satisfaying the above conditions
    % we have to add an upstream arc (speed w) and an upstream node
    % warning : this upstream arc is downstream on the arterial
    NumNodes = NumNodes + 1; % new node
    NumArcs = NumArcs + 1; % new arc
    % --- (a) we complete the information about the current node ---
    Node(Liste(i_node)).UpstreamArc(end+1) = NumArcs;
    Node(Liste(i_node)).NumUpstreamArc = Node(Liste(i_node)).NumUpstreamArc + 1;
    CurrentNode = Node(Liste(i_node));
    % --- (b) we add the final node ---
    Node(NumNodes).T = CurrentNode.T - (Link(end).Position - CurrentNode.X) / FD.w; % time corresponding to the next point
    Node(NumNodes).X = Link(end).Position; % position of the new node : start of the arterial
    Node(NumNodes).UpstreamArc = [];
    Node(NumNodes).NumUpstreamArc = 0;
    Node(NumNodes).DownstreamArc = NumArcs; % id of the new arc
    Node(NumNodes).NumDownstreamArc = 1;
    Node(NumNodes).RankX = RankXMax + 1; % convention : end of the arterial
    Node(NumNodes).Type = 1 ; % origin of a path (necessary for begin a path)
    NewNode = Node(NumNodes);
    % --- (c) we add the final arc ---
    Arc(NumArcs).Pamont = [NewNode.T NewNode.X];
    Arc(NumArcs).Paval = [CurrentNode.T CurrentNode.X];
    Arc(NumArcs).RankX = 0; % Convention: an oblique arc has a rankx equal to 0
    Arc(NumArcs).CostRate = Link(RankXMax + 1).NumLanes * FD.w * FD.kappa; % wave at speed w
    Arc(NumArcs).duration=[Arc(NumArcs).Paval(1) - Arc(NumArcs).Pamont(1)];
    Arc(NumArcs).UpstreamNode = NumNodes;
    Arc(NumArcs).DownstreamNode = Liste(i_node);
end

% ----- (4) Completion of descending paths on the first link -----
Condition1 = ([Node.RankX] == 1); % nodes on the first signal
Condition2a = ([Node.Type] == 1); % ends of red phases
Condition2b = ([Node.Type] == 3); % ends of arcs at the speed w
Liste = find(Condition1 & (Condition2a | Condition2b));
for i_node = 1 : length(Liste)
    % for each node satisfaying the above conditions
    % we have to add a downstream arc (speed w) and a downstream node
    % warning : this downstream arc is upstream on the arterial
    NumNodes = NumNodes + 1; % new node
    NumArcs = NumArcs + 1; % new arc
    % --- (a) we complete the information about the current node ---
    Node(Liste(i_node)).DownstreamArc(end+1) = NumArcs;
    Node(Liste(i_node)).NumDownstreamArc = Node(Liste(i_node)).NumDownstreamArc + 1;
    CurrentNode = Node(Liste(i_node));
    % --- (b) we add the final node ---
    Node(NumNodes).T = CurrentNode.T - (0 - CurrentNode.X) / FD.w; % time corresponding to the next point
    Node(NumNodes).X = 0; % position of the new node : start of the arterial
    Node(NumNodes).UpstreamArc = NumArcs; % id of the new arc
    Node(NumNodes).NumUpstreamArc = 1;
    Node(NumNodes).DownstreamArc = [];
    Node(NumNodes).NumDownstreamArc = 0;
    Node(NumNodes).RankX = 0; % convention : start of the arterial
    Node(NumNodes).Type = 3 ; % end of a wave at speed w
    NewNode = Node(NumNodes);
    % --- (c) we add the final arc ---
    Arc(NumArcs).Pamont = [CurrentNode.T CurrentNode.X];
    Arc(NumArcs).Paval = [NewNode.T NewNode.X];
    Arc(NumArcs).RankX = 0; % Convention: an oblique arc has a rankx equal to 0
    Arc(NumArcs).CostRate = Link(1).NumLanes * FD.w * FD.kappa; % wave at speed w
    Arc(NumArcs).duration=[Arc(NumArcs).Paval(1) - Arc(NumArcs).Pamont(1)];
    Arc(NumArcs).UpstreamNode = Liste(i_node);
    Arc(NumArcs).DownstreamNode = NumNodes;
end
end

function [ Cut ] = shortest_variational_paths( Link, Arc, Node, CalculationMethod, G )

% Link structure : defined the link main characteristics
    % Link.ID : Link ID
    % Link.Length : Link length [m]
    % Link.Position : Global position of traffic signal [m] (auto)
    % Link.NumLanes : Number of Lanes
    % Link.SignalGT : Green Time [s]
    % Link.SignalC  : Signal cycle [s]
    % Link.SignalDelay : lag for the signal becoming green inside the cycle [s]
    % Link.ListRedTime : List of times when the traffic signal turns red [s] (auto)
    % Link.ListGreenTime : List of times when the traffic signal turns green [s] (auto)
    % Link.Capacity : Downstream capacity due to traffic signal [veh/s] (auto)
% Arc structure : Defined the arc main characteristics
    % Arc.Pamont : Coordinates of the upstream point of the arc
    % Arc.Paval : Coordianates of the downstream point of the arc
    % Arc.RankX : Associated link ID for horizontal arcs
    % Arc.CostRate : CostRate associated to the arc (The total cost is equal to the costrate multiplied by the time distance of the arc
    % Arc.UpstreamNode : Node ID of the upstream point
    % Arc.DownstreamNode : Node ID of the downstream point
% Node structure: Defined the node main characteristics and connections to arcs
    % Node.T : T-Coordinate of the nodes (t,x) plane
    % Node.X : X-oordinate of the nodes (t,x) plane
    % Node.UpstreamArc : List of Id of the upstream arcs
    % Node.NumUpstreamArc : Number of upstream arcs
    % Node.DownstreamArc : List of Id of the downstream arcs
    % Node.NumDownstreamArc : Number of donwstream arcs
    % Node.Type:
        % 0, the node is connected with a downstream red phase;
        % 1, the node is connected with a downstream green phase;
        % 2 the node is a intermediate node created by intersecting oblique arc and wille be the origin of a wave at speed u
        % 3 the node is a intermediate node created by intersecting oblique arc and wille be the origin of a wave at speed w
% Cut structure: Defined the cut main characteristics and path
    % Cut.InitialNode : initial node of the path
    % Cut.FinalNode : final node of the path
    % Cut.Speed : observater mean speed
    % Cut.r : overtaking flow on this path
    % Cut.TravelTime : travel time of the path
    % Cut.Path : Vector containing the id of the traveled nodes

% ----- Defintion of the CostMatrix for the shortest path problem -----
NumLinks = length(Link);
NumNodes = length(Node);
NumArcs = length(Arc);

% Temp_AdjMatrix is a NxN adjacency matrix, where A(I,J) is nonzero (=1)
% if and only if an edge connects point I to point J
Temp_AdjMatrix = zeros(NumNodes,NumNodes); % default value if two nodes are not linked
% Temp_CostMatrix is a NxN cost (perhaps distance) matrix, where C(I,J)
% contains the value of the cost to move from point I to point J
Temp_CostMatrix = +Inf * ones(NumNodes,NumNodes); % default value if two nodes are not linked

for i=1:NumArcs
    % for each couple of nodes linked together
    Temp_AdjMatrix(Arc(i).UpstreamNode,Arc(i).DownstreamNode) = 1; % these both nodes are linked together
    Temp_CostMatrix(Arc(i).UpstreamNode,Arc(i).DownstreamNode) = Arc(i).CostRate * (Arc(i).Paval(1) - Arc(i).Pamont(1)); % cost = unitary cost * arc duration
end

% ----- Solving the shortest path problem -----
Cut = struct;
IdCut=0; % initialisation

% ---- for positive mean speeds of observers u
% --- nodes to investigate
if CalculationMethod.BoundaryPaths
    ListInitialNodeId = find( ([Node(:).RankX]==0) & ([Node(:).T]>=0) & ([Node(:).Type]==1)); % nodes at the most upstream signal, at end of a red phase, and after the "start" of all other signals
    if CalculationMethod.EndRed
        FinalNodeId=find(([Node(:).RankX]==NumLinks) & ([Node(:).Type]==1)); % nodes at the last signal (most downstream). We can want it to be at end of red phase
    else
        FinalNodeId=find([Node(:).RankX]==NumLinks); % nodes at the last signal (most downstream)
    end
else
    ListInitialNodeId = find( ([Node(:).RankX]==1) & ([Node(:).T]>=0) & ([Node(:).Type]==1)); % nodes at the most upstream signal, at end of a red phase, and after the "start" of all other signals
    if CalculationMethod.EndRed
        FinalNodeId=find(([Node(:).RankX]==NumLinks-1) & ([Node(:).Type]==1)); % nodes at the last signal (most downstream). We can want it to be at end of red phase
    else
        FinalNodeId=find([Node(:).RankX]==NumLinks-1); % nodes at the last signal (most downstream).
    end
end

for i_initial = 1 : length(ListInitialNodeId)
    InitialNodeId = ListInitialNodeId(i_initial); % considered initial node
    % --- find the shortest-paths
    shortestDists = distances(G, InitialNodeId,'Method','acyclic');
    totalcost = shortestDists(FinalNodeId);
    
    % --- compute the existing cuts
    for i=1:1:length(FinalNodeId)
        % for all investigated final nodes
        Temp_TravelTime = Node(FinalNodeId(i)).T - Node(InitialNodeId).T;
        if totalcost(i) < +Inf
            % if the observer can reach the investigated node, there is a cut
            IdCut = IdCut+1; % new cut
            Cut(IdCut).InitialNode = InitialNodeId; % initial node
            Cut(IdCut).FinalNode = FinalNodeId(i); % current investigated node
            Cut(IdCut).Speed = (Node(FinalNodeId(i)).X - Node(InitialNodeId).X) / Temp_TravelTime; % observater mean speed
            Cut(IdCut).r = totalcost(i) / Temp_TravelTime; % overtaking flow
            Cut(IdCut).TravelTime = Temp_TravelTime;
        end
    end
end

% [P,d] = shortestpath(G,254,281)


% ---- for negative mean speeds of observers
% --- nodes to investigate
if CalculationMethod.BoundaryPaths
    ListInitialNodeId = find( ([Node(:).RankX]==NumLinks) & ([Node(:).T]>=0) & ([Node(:).Type]==1));  % node at the most downstream signal, at end of a red phase, and after the "start" of all other signals
    if CalculationMethod.EndRed
        FinalNodeId=find(([Node(:).RankX]==0) & ([Node(:).Type]==1)); % nodes at the first signal (most upstream). We can want it to be at end of red phase
    else
        FinalNodeId=find([Node(:).RankX]==0); % nodes at the first signal (most upstream).
    end
else
    ListInitialNodeId = find( ([Node(:).RankX]==NumLinks-1) & ([Node(:).T]>=0) & ([Node(:).Type]==1));  % node at the most downstream signal, at end of a red phase, and after the "start" of all other signals
    if CalculationMethod.EndRed
        FinalNodeId=find(([Node(:).RankX]==1) & ([Node(:).Type]==1)); % nodes at the first signal (most upstream). We can want it to be at end of red phase
    else
        FinalNodeId=find([Node(:).RankX]==1); % nodes at the first signal (most upstream).
    end
end

for i_initial = 1 : length(ListInitialNodeId)
    InitialNodeId = ListInitialNodeId(i_initial); % considered initial node
    % --- find the shortest-paths
    shortestDists = distances(G, InitialNodeId,'Method','acyclic');
    totalcost = shortestDists(FinalNodeId);

    % --- compute the existing cuts
    for i=1:1:length(FinalNodeId)
        % for all investigated final nodes
        Temp_TravelTime = Node(FinalNodeId(i)).T - Node(InitialNodeId).T;
        if totalcost(i) < +Inf
            % if the observer can reach the investigated node, there is a cut
            IdCut = IdCut+1; % new cut
            Cut(IdCut).InitialNode = InitialNodeId; % initial node
            Cut(IdCut).FinalNode = FinalNodeId(i); % current investigated node
            Cut(IdCut).Speed = (Node(FinalNodeId(i)).X - Node(InitialNodeId).X) / Temp_TravelTime; % observater mean speed
            Cut(IdCut).r = totalcost(i) / Temp_TravelTime; % overtaking flow
            Cut(IdCut).TravelTime = Temp_TravelTime;
        end
    end
end
    

% ---- for null mean speeds of observers
NumLinks = length(Link) - 1;
for i_link = 1 : NumLinks
    % --- nodes to investigate
    Tini = Link(i_link).ListGreenTime(find( Link(i_link).ListGreenTime>=0,1,'first'));
    Tend = Link(i_link).ListGreenTime(find( Link(i_link).ListGreenTime<=CalculationMethod.TimeWindow,1,'last'));
    InitialNodeId = find( ([Node(:).RankX]==i_link) & ([Node(:).T]==Tini));  % first node at the most downstream signal, at end of a red phase, and after the "start" of all other signals
    FinalNodeId = find( ([Node(:).RankX]==i_link) & ([Node(:).T]==Tend)); % node at the first signal (most upstream). We can want it to be at end of red phase
    %InitialNodeId = find( ([Node(:).RankX]==i_link) & ([Node(:).T]>=0) & ([Node(:).Type]==1),1,'first');  % first node at the most downstream signal, at end of a red phase, and after the "start" of all other signals
    %FinalNodeId=find([Node(:).RankX]==i_link & [Node(:).T]>Node(InitialNodeId).T & [Node(:).Type]==1,1,'last'); % node at the first signal (most upstream). We can want it to be at end of red phase
    
    % --- find the shortest-paths
    shortestDists = distances(G, InitialNodeId,'Method','acyclic');
    totalcost = shortestDists(FinalNodeId);

    % --- compute the existing cuts
    Temp_TravelTime = Node(FinalNodeId).T - Node(InitialNodeId).T;
    if totalcost < +Inf
        % if the observer can reach the investigated node, there is a cut
        IdCut = IdCut+1; % new cut
        Cut(IdCut).InitialNode = InitialNodeId; % initial node
        Cut(IdCut).FinalNode = FinalNodeId; % current investigated node
        Cut(IdCut).Speed = 0; % observer mean speed
        Cut(IdCut).r = totalcost / Temp_TravelTime; % overtaking flow
        Cut(IdCut).TravelTime = Temp_TravelTime;
    end
end
end

function [ MacroFundDiag, ListLowerCuts ] = mfd_by_lower_cuts( Cut )
% Selection of the lowerCut that defines the MFD

% Cut structure: Define the cut main characteristics and path
    % Cut.InitialNode : initial node of the path
    % Cut.FinalNode : final node of the path
    % Cut.Speed : observater mean speed
    % Cut.r : overtaking flow on this path
    % Cut.TravelTime : travel time of the path
    % Cut.Path : Vector containing the id of the traveled nodes
% MacroFundDiag structure: Define the characteristic points of the MFD
    % MacroFundDiag.k : vector of concentration
    % MacroFundDiag.q : vector of flow
    % MacroFundDiag.LeftCut : vector with the id of cuts on the left side of the point
    % MacroFundDiag.RightCut : vector with the id of cuts on the right side of the point
% ListLowerCuts is a vector containing the id of lower cuts

% ----- Cuts to consider -----
% Creation of VRaS : containing the overtaking flow and the speed of all cuts
VRaS(:,1) = [Cut.r]'; % overtaking flow on each path
VRaS(:,2) = [Cut.Speed]'; % observer mean speed for each path
VRaS(:,3) = (1:length(Cut))'; % Cut id

% Check for standard deviation mean ratio of costs for cuts.
VRaS2(:,1) = [Cut.Speed]'; % overtaking flow on each path
VRaS2(:,2) = [Cut.r]'; % observer mean speed for each path
VRaS2(:,3) = (1:length(Cut))'; % Cut id
VRaS2 = sortrows(VRaS2);
stdCuts = calcStandardErrorCuts(VRaS2);

VRaS = sortrows(VRaS); % sort of the table from the smallest to highest overtaking flows, and from the smallest to the highest speed

%%%%%
% HERE I MAKE SURE THE MEAN STATIONARY CUT IS TAKEN.
test = VRaS(VRaS(:,2)==0,1);
VRaS(VRaS(:,2)==0,1)=mean(test);
%%%%%

[~,list,~] = unique(VRaS(:,1),'first'); % id of overtaking flows with the smallest speed (tight cut)
VRaS = VRaS(list,:); % we keep only different overtaking flows with the smallest speed (tight cut)
NumCuts = size(VRaS,1);

% ----- Creation of MFD by intersections between cuts -----
KsuchQis0 = -VRaS(:,1)./VRaS(:,2); % intersections between the cuts and the line Q = 0
Kmax = min(KsuchQis0(KsuchQis0>0)); % tight cut for congested part at Q = 0
Qmin = min(VRaS(:,1)); % tight cut for uncongested part at K = 0
FirstLoweverCutId = find(VRaS(:,1) == Qmin,1,'first'); % id of the tight cut for congested part at Q = 0
LastLoweverCutId = find(KsuchQis0==min(KsuchQis0(KsuchQis0>0)) & KsuchQis0>0,1,'first'); % id of the tight cut for congested part at Q = 0
% the both extremal cuts are the first one in VRaS, and LastLoweverCutId in VRaS
% Their intersection point is:
FirstIntersection.k = (VRaS(LastLoweverCutId,1) - VRaS(1,1)) / (VRaS(1,2) - VRaS(LastLoweverCutId,2));
FirstIntersection.q = VRaS(1,1) + VRaS(1,2) * FirstIntersection.k;

% --- initialisation of the algorithm to find the lower cuts
MacroFundDiag.k = [0; FirstIntersection.k; Kmax]; % vector of density
MacroFundDiag.q = [Qmin; FirstIntersection.q; 0]; % vector of flow
MacroFundDiag.LeftCut = [0; FirstLoweverCutId; LastLoweverCutId]; % vector with the id of cuts on the left side of the point
MacroFundDiag.RightCut = [FirstLoweverCutId; LastLoweverCutId; 0]; % vector with the id of cuts on the right side of the point

% --- Algorithm to find the lower cuts
for i_cut = 2 : NumCuts
    % for each cut
    % equation of this cut: y = VRaS(i_cut,1) + VRaS(i_cut,2) . x
    Q_of_Cut_i = VRaS(i_cut,1) + VRaS(i_cut,2) * MacroFundDiag.k; % ordinate of the cut at temporary points of the MFD
    list = find( (Q_of_Cut_i - MacroFundDiag.q) < -10^-5); % list of point that won't be more MFD points
    if ~isempty(list)
        % so this cut is a lower cut
        LeftCut = MacroFundDiag.LeftCut(list(1)); % left-cut intersected by the new cut
        LeftIntersection = intersection( VRaS(LeftCut,2), VRaS(LeftCut,1), VRaS(i_cut,2), VRaS(i_cut,1) ); % intersection of this cut with the temporary MFD on the left
        RightCut = MacroFundDiag.RightCut(list(end)); % right-cut intersected by the new cut

        RightIntersection = intersection( VRaS(RightCut,2), VRaS(RightCut,1), VRaS(i_cut,2), VRaS(i_cut,1) );  % intersection of this cut with the temporary MFD on the right        
        % MFD changes
        MacroFundDiag.k = [MacroFundDiag.k(1:list(1)-1); LeftIntersection.k; RightIntersection.k; MacroFundDiag.k(list(end)+1:end)];
        MacroFundDiag.q = [MacroFundDiag.q(1:list(1)-1); LeftIntersection.q; RightIntersection.q; MacroFundDiag.q(list(end)+1:end)];
        MacroFundDiag.LeftCut = [MacroFundDiag.LeftCut(1:list(1)-1); LeftCut; i_cut; MacroFundDiag.LeftCut(list(end)+1:end)];
        MacroFundDiag.RightCut = [MacroFundDiag.RightCut(1:list(1)-1); i_cut; RightCut; MacroFundDiag.RightCut(list(end)+1:end)];
    end % else, the MFD is not changed by this cut
end

% --- we ensure that all points are not the same
Tmp_pts_length = length(MacroFundDiag.k);
i = 1;
while i < Tmp_pts_length
    if abs(MacroFundDiag.k(i+1) - MacroFundDiag.k(i)) < 10^-5
        % these both points are the same
        MacroFundDiag.k = MacroFundDiag.k([1:i,i+2:end]);
        MacroFundDiag.q = MacroFundDiag.q([1:i,i+2:end]);
        MacroFundDiag.LeftCut = MacroFundDiag.LeftCut([1:i,i+2:end]);
        MacroFundDiag.RightCut = MacroFundDiag.RightCut([1:i-1,i+1:end]);
        % the second point is deleted
        Tmp_pts_length = Tmp_pts_length - 1;
    else
        i = i + 1;
    end
end

% --- we keep the lower cuts id
ListLowerCuts = VRaS(MacroFundDiag.LeftCut(2:end),3); % id of lower cuts
MacroFundDiag.LeftCut(2:end) = ListLowerCuts; % id of lower cuts
MacroFundDiag.RightCut(1:end-1) = ListLowerCuts; % id of lower cuts

end

function [ Int ] = intersection( a,b,c,d )
% function to compute the interesection between two lines of equations:
    % y = ax + b
    % y = cx + d

if a~=c
    % different speeds
    Int.k = (d-b) / (a-c);
    Int.q = a*Int.k + b;
else
    % same speeds
    disp('error in ''intersection'' in ''mfd_by_lower_cuts'' !')
    disp([a,b,c,d]);
end

end

function [temp2] = calcStandardErrorCuts(VRaS2)

% round speeds to first nachkommastelle
temp = [round(VRaS2(:,1),1) VRaS2(:,2)];

% get unique speeds
temp2 = unique(temp(:,1));

% allocate space for mean/std/ratio
temp2 = [temp2 zeros(length(temp2),3)];

% calc mean/std/ratio for each unique speed
for i=1:1:length(temp2)
    temp2(i,2) = mean(temp(temp(:,1) == temp2(i),2));
    temp2(i,3) = std(temp(temp(:,1) == temp2(i),2));
    temp2(i,4) = temp2(i,3)/temp2(i,2)*100;
end


end

function [G] = plot_graph(Arc, Node)
s = [Arc.UpstreamNode];
t = [Arc.DownstreamNode];
weights = [Arc.CostRate].*[Arc.duration];
% weights = [Arc.CostRate];
G = digraph(s,t,weights);
% p = plot(G,'XData',[Node.T],'YData',[Node.X],'EdgeLabel',G.Edges.Weight)
end