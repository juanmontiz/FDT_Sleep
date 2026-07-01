
%%% A2_SIGMA_EXPLORATION_ANALYZE_RESULTS  Plots D(sigma) curves; identifies optimal sigma.
% Author: Juan Manuel Monti
%
% Loads sigmaexplore_Sleep_COND_*.mat and plots mean quadratic error D(sigma)
% with shaded std for Wakefulness (blue) and Deep Sleep N3 (red). Marks the
% sigma minimizing D with a dashed vertical line and a marker dot.
%
% Generates three figure sets:
%   sigma_variartion_Sleep_COND_1_NSIM_<N>.{svg,png}      (W only)
%   sigma_variartion_Sleep_COND_2_NSIM_<N>.{svg,png}      (N3 only)
%   sigma_variartion_Sleep_CONDS_1-2_NSIM_<N>.{svg,png}   (W and N3 overlaid)
%
% Requires: sigmaexplore_Sleep_COND_*.mat, shadedErrorBar

clearvars

NSIM = 100

%%%% TS NORMALIZATION
% z-score std --> zcs_std

%%%% DIFFERENCE BETWEEN C(t,t') MATRICSE
% norm --> norm(Cemp - Csim)
% quad --> mean(mean((Cemp - Csim)^2))
%  d1 -->  mean(mean(abs(Cemp - Csim)))
%  d2 -->  sqrt(mean(mean((Cemp - Csim)^2))


%%% COND 1 %%%
load(sprintf('sigmaexplore_Sleep_COND_1_NSUB_15_NSUBSIM_%d_DFILT_0_FREQSUB_0.mat',NSIM))
%%% z-score
%  e1 = errorC_norm_sim_zsc_std;
e1 = errorC_quad_sim_zsc_std;
%  e1 = errorC_d1_sim_zsc_std;
%  e1 = errorC_d2_sim_zsc_std;
sigma1 = sigma;

%%% COND 2 %%%
load(sprintf('sigmaexplore_Sleep_COND_2_NSUB_15_NSUBSIM_%d_DFILT_0_FREQSUB_0.mat',NSIM))
%%% z-score
%  e2 = errorC_norm_sim_zsc_std;
e2 = errorC_quad_sim_zsc_std;
%  e2 = errorC_d1_sim_zsc_std;
%  e2 = errorC_d2_sim_zsc_std;
sigma2 = sigma;

clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COND 1 ONLY
minx = 0;
maxx = 0.5;
% miny = 0.965;
% maxy = .99;
set(gca,'Xtick',0:0.05:1.5)
xlim([minx maxx])
% ylim([miny maxy])
% xlabel('\sigma'),ylabel('mean(Cemp - Csim)^2')
xlabel('\sigma'),ylabel('D')
box on

aux = e1;
sigma = sigma1(1:length(aux(:,1)));
col = 'b';
dispaux1 = 'W';
%% Shadow
stdaux = std(aux,[],2);
shadedErrorBar(sigma,mean(aux,2),stdaux,'lineprops',col),hold on
% vertical line
sigmaux1 = sigma(find(mean(aux,2)==min(mean(aux,2))));
xline(sigmaux1,'--','color',col,'LineWidth',1,'DisplayName',num2str(sigmaux1))
% dot
plot(sigmaux1,min(mean(aux,2)),'o','MarkerEdgeColor','k','MarkerFaceColor', col)

legend(dispaux1,['\sigma = ',num2str(sigmaux1)],...
       'Location','NorthEast')
hold off

%%% save figure
f = gcf;
file = sprintf('sigma_variartion_Sleep_COND_1_NSIM_%d.svg',NSIM);
print(f,'-dsvg',file)
file = sprintf('sigma_variartion_Sleep_COND_1_NSIM_%d.png',NSIM);
saveas(f,file)
close(f)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COND 2 ONLY
clf
minx = 0;
maxx = 0.5;
% miny = 0.965;
% maxy = .99;
set(gca,'Xtick',0:0.05:1.5)
xlim([minx maxx])
% ylim([miny maxy])
% xlabel('\sigma'),ylabel('mean(Cemp - Csim)^2')
xlabel('\sigma'),ylabel('D')
box on

aux = e2;
sigma = sigma2(1:length(aux(:,1)));
col = 'r';
dispaux2 = 'N3';
%% Shadow
stdaux = std(aux,[],2);
shadedErrorBar(sigma,mean(aux,2),stdaux,'lineprops',col),hold on
% vertical line
sigmaux2 = sigma(find(mean(aux,2)==min(mean(aux,2))));
xline(sigmaux2,'--','color',col,'LineWidth',1,'DisplayName',num2str(sigmaux2))
% dot
plot(sigmaux2,min(mean(aux,2)),'o','MarkerEdgeColor','k','MarkerFaceColor', col)

legend(dispaux2,['\sigma = ',num2str(sigmaux2)], ...
       'Location','NorthEast')
hold off

%%% save figure
f = gcf;
file = sprintf('sigma_variartion_Sleep_COND_2_NSIM_%d.svg',NSIM);
print(f,'-dsvg',file)
file = sprintf('sigma_variartion_Sleep_COND_2_NSIM_%d.png',NSIM);
saveas(f,file)

close(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
%%% COND 1 & 2
%% COND 1
minx = 0;
maxx = 0.5;
%  miny = 0.92;
%  maxy = 1;
set(gca,'Xtick',0:0.05:1.5)
xlim([minx maxx])
%  ylim([miny maxy])
%  xlabel('\sigma'),ylabel('mean(Cemp - Csim)^2')
xlabel('\sigma'),ylabel('D')
box on

sigma = sigma1;
aux = e1;
col = 'b';
dispaux1 = 'W';
%% Shadow
stdaux = std(aux,[],2);
shadedErrorBar(sigma,mean(aux,2),stdaux,'lineprops',col),hold on
% vertical line
sigmaux1 = sigma(find(mean(aux,2)==min(mean(aux,2))));
xline(sigmaux1,'--','color',col,'LineWidth',1,'DisplayName',num2str(sigmaux1))
% dot
plot(sigmaux1,min(mean(aux,2)),'o','MarkerEdgeColor','k','MarkerFaceColor', col)

%% COND 2
aux = e2;
sigma = sigma2(1:length(aux(:,1)));
col = 'r';
dispaux2 = 'N3';
%% Shadow
stdaux = std(aux,[],2);
shadedErrorBar(sigma,mean(aux,2),stdaux,'lineprops',col),hold on
% vertical line
sigmaux2 = sigma(find(mean(aux,2)==min(mean(aux,2))));
xline(sigmaux2,'--','color',col,'LineWidth',1,'DisplayName',num2str(sigmaux2))
% dot
plot(sigmaux2,min(mean(aux,2)),'o','MarkerEdgeColor','k','MarkerFaceColor', col)

legend(dispaux1,['\sigma = ',num2str(sigmaux1)],'', ...
       dispaux2,['\sigma = ',num2str(sigmaux2)],'', ...
       'Location','NorthEast')
hold off

%%% save figure
f = gcf;
file = sprintf('sigma_variartion_Sleep_CONDS_1-2_NSIM_%d.svg',NSIM);
print(f,'-dsvg',file)
file = sprintf('sigma_variartion_Sleep_CONDS_1-2_NSIM_%d.png',NSIM);
saveas(f,file)


close(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
