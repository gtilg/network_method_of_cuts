function [mfd] = estimateCapacityBranch(mfd, network)
%ESTIMATECAPACITYBRANCH Estimates the capacity branch of the MFD
%   Detailed explanation goes here
q = mfd.q;
hyperlinks = network.links(:,{'id','length','cycle','green','maxflow_i'});
hyperlinks = unique(hyperlinks,'rows','first');

hyperlinks.capacity = [hyperlinks.maxflow_i] .* [hyperlinks.green]./[hyperlinks.cycle].*[hyperlinks.length];
capacity_net = sum([hyperlinks.capacity])/sum([hyperlinks.length]); % length-weighted average

q(q>capacity_net*1800)=capacity_net*1800;
mfd.q = q;

end

