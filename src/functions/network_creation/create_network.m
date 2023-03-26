function [network] = create_network(tr_scenario)

network = struct();

rng(tr_scenario)
% rng(0);
rand(1,1);

load('links_sf_bi.mat');
load('conn_sf_bi.mat');

% Deviate turning ratios a bit
connections = sortrows(connections,'fromLink','ascend');
for i=1:2:height(connections)
    connections.turningRatio(i) = round(0.25 + rand(1,1)*0.5,2);
    connections.turningRatio(i+1) = 1-connections.turningRatio(i);
end

% Sort links according to direction of travel
links = sortrows(links,{'corridor','order'},'ascend');

nLinks = height(links);

% find connections according to turning ratios
for i=1:nLinks
    cands = connections(links.index(i) == connections.toLink,:);

    if height(cands)==2
        links.j_link(links.index == cands.fromLink(cands.fromCorr==cands.toCorr)) = cands.fromLink(cands.fromCorr~=cands.toCorr);
    end
    
    for j=1:height(cands)
        if links.corridor(links.index==cands.fromLink(j)) == links.corridor(links.index==cands.toLink(j))
            links.alpha_ij(links.index==cands.fromLink(cands.fromLink == cands.fromLink(j))) = 1-cands.turningRatio(j); % find alpha_ij
            
        elseif links.corridor(links.index==cands.fromLink(j)) ~= links.corridor(links.index==cands.toLink(j))
            if ~isempty(cands.fromLink(cands.fromLink ~= cands.fromLink(j))) % when no turning ratios are specified in input file
                links.alpha_ji(links.index==cands.fromLink(cands.fromLink ~= cands.fromLink(j))) = cands.turningRatio(j); % find alpha_ji
            end
        end
        
    end
end

% add origin and destination
links.origin = zeros(nLinks,1);
links.destination = zeros(nLinks,1);
for i=1:nLinks
    if links.fromIntersection(i) > 1000
        links.origin(i) = 1;
    end
    if links.toIntersection(i) > 1000
        links.destination(i) = 1;
    end
end

% Find upstream i and j
for i=1:nLinks
    if links.origin(i)~=1
        links.upstream_i(i) = links.index(i-1);
        links.upstream_j(i) = links.j_link(i-1);
    else
        links.upstream_i(i) = NaN;
        links.upstream_j(i) = NaN;
    end
end

% Find downstream i and j
for i=1:nLinks
    if links.destination(i)~=1
        links.downstream_i(i) = links.index(i+1);
        links.downstream_j(i) = links.j_link(i+1);
    else
        links.downstream_i(i) = NaN;
        links.downstream_j(i) = NaN;
    end
end

% add maxflows
links.maxflow_i = ones(nLinks,1);
links.maxflow_j = ones(nLinks,1);

% add cycle, green and offset
links.green = ones(nLinks,1)*45;
links.cycle = ones(nLinks,1)*90;
links.red = links.cycle - links.green;
links.offset = zeros(nLinks,1);
links.offset(links.corridor > 5 & links.corridor < 13)=45;
links.offset(links.corridor > 17)=45;

% Delete unnecessary variables
columnIndicesToDelete = [6];
links(:,columnIndicesToDelete)=[];

% Rename if necessary
links = renamevars(links,["index"], ["id"]);

% Reordering
links = links(:,{'id' 'length' 'cumLength' 'corridor' 'cycle' 'green' 'offset' 'red' 'alpha_ij' 'alpha_ji' 'j_link' 'origin' 'destination' 'upstream_i' 'upstream_j' 'downstream_i' 'downstream_j'});

network.links = links;

%% Connections
% Rename if necessary
connections = renamevars(connections,["turningRatio", "index"], ["Share", "id"]);

% Reordering
connections = connections(:,{'id' 'intersection' 'fromLink' 'toLink' 'Share' 'fromCorr' 'toCorr'});

network.connections = connections;

end
