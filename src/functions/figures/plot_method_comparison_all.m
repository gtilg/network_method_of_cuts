function plot_method_comparison_all(mfd_net_nvt, mfd_net_approx1, mfd_net_approx2)
%PLOT_FIGURES Plots the MFDs from all evaluated methods
%   This refers to the CTM ground truth data, the original method of cuts
%   by Daganzo & Geroliminis (DG), by Leclercq and Geroliminis (LG), and by
%   Aghamohammadi and Laval (AL).

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Load data
% CTM: Sioux Falls
load('ctm_results_7200_1.mat')
k_groundtruth_ctm = mfd_k;
q_groundtruth_ctm = mfd_q;

% Original DG data
load('origMC_DG.mat')
k_daganzo = k;
q_daganzo = q;

% Original LG data
load('origMC_LG.mat')
k_lec = mfd_lec.k;
q_lec = mfd_lec.q;

% Original AL data
load('origMC_AL.mat')
k_lav = mfd_lav.k;
q_lav = mfd_lav.q;

%% Plot nMC/MC MFDs
figure;
scatter(k_groundtruth_ctm,q_groundtruth_ctm,'d')
hold on
plot(k_daganzo, q_daganzo,'-', 'Color', [0.7 0.7 0.7])
hold on
plot(k_lec, q_lec,'--', 'Color', [0.7 0.7 0.7])
hold on
plot(k_lav, q_lav,':', 'Color', [0.7 0.7 0.7])
hold on
plot(mfd_net_nvt.k,mfd_net_nvt.q,'-k')
hold on
plot(mfd_net_approx2.k,mfd_net_approx2.q,'--k')
hold on
plot(mfd_net_approx1.k,mfd_net_approx1.q,':k')
box on

% Polish the figure
ylim([0 1000])
xlim([0 175])
xticks([0 25 50 75 100 125 150 175])
yticks([0 200 400 600 800 1000])
xlabel('K [veh/km]')
ylabel('Q [veh/h]')
legend('CTM', 'DG', 'LG', 'AL', 'nVT', 'FS-QP', 'LS-QP', 'Location','northeast','Color','none','Orientation','vertical')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gca,'linewidth',1.2)
set(gcf,'Position', [10 10 800 800])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gca, 'Type', 'Scatter'),'SizeData', 30,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);

width = 6.5; % inches
nrow = 1;
height = 6.5; % inches
set(gcf,'PaperUnits','inches')
set(gcf,'PaperSize',[width height],'PaperPosition',[0 0 width height])
print('-dpdf','figures/result.pdf')
end

