function [mfd] = estimateCongestedBranch(mfd, kappa_net, FD)
%ESTIMATECONGESTEDBRANCH Estimates the congetsed branch of the MFD
%   Relying on the canonical transformation of the FD (Laval, Castrillon;
%   2015), we transform the free-flow branch and then mirror the data at
%   kappa_net/2.

q = mfd.q;
k = mfd.k;

% Find start of backward wave cut
q_diff = diff(q);
id_bw_start = find(q_diff<0,1)+1;
k(id_bw_start:end)=[];
q(id_bw_start:end)=[];

% Transform with kappa_net.
w_net = FD.qmax/(kappa_net - FD.kc);
% w_net = FD.w;
k_prime = k/kappa_net/1000 - 1/2 * (1 - (1/w_net - 1/FD.u)*q/FD.qmax/3600);
q_new = q(k_prime<=0);
k_prime(k_prime>0)=[];

k_prime = [k_prime, -flip(k_prime)];
q_prime = [q_new,flip(q_new)];

for i = 1:length(k_prime)
    k_trans(i)=(k_prime(i)+0.5*(1-(1/w_net-1/FD.u)*q_prime(i)/FD.qmax/3600))*kappa_net*1000;
end

mfd.q = [0 q_prime 0];
mfd.k = [0 k_trans kappa_net*1000];

end

