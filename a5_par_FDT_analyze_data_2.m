%%% A5_PAR_FDT_ANALYZE_DATA_2  Plots parcel-wise FDT metrics: Wakefulness vs Deep Sleep.
% Author: Juan Manuel Monti
%
% Loads metrics_PAR_*_trapz.mat for both conditions and generates violin plots
% and line plots comparing W vs N3 across parcels for xMi, xMd, Mi, Md, MX, xMX.
% Also exports CSV files for Python-based violin plots (a6_boxplots_M_v2.ipynb).
%
% Figures saved per metric:
%   Sleep_NSIM_<N>_<metric>_viol.{svg,png}  violin plot (W vs N3)
%   Sleep_NSIM_<N>_<metric>_plot.{svg,png}  per-parcel line plot
%
% CSV exports (always active):
%   Sleep_NSIM_<N>_<metric>_par_for_violin.csv  columns: cond, value
%
% Requires: metrics_PAR_Sleep_COND_*_trapz.mat, violinplot

clearvars

%  Nsimulations = [100 250 750 1000 2500 7500 10000]
Nsimulations = [10000]
for NSUBSIM = Nsimulations
%%% Number of simulations for each subject
NSUBSIM

%%% Choose dataset and other parameters %%%
%%% SLEEP DATASET

%%% SET CONDITION
% 1 --> Awake
% 2 --> Sleep

%%% Number of subjects
NSUB = 15;

%%% Use filtered or non-filtered data (ts, force, noise) to calculate C, A and R
%   DATAFILTER = 1 --> use filtered data
%   DATAFILTER = 0 --> use non-filtered data
DATAFILTER = 0;

%%% Linearization: Linear (LIN=1 NOT YET IMPLEMENTED) or Cubic (LIN=0) terms in Hopf
LIN = 0;

%%% Hopf initialization
%   HOPFINT = 1 --> initializes for each individual simulation
%   HOPFINT = 0 --> initializes ONLY for the first simulation (of each subject)
%                   and all subsequent use the same z0
HOPFINIT = 0;

%%% Hopf Frequencies
%   FREQSUB = 1 --> different frequencies for each subject
%   FREQSUB = 0 --> frequencies for the mean Power Spectra
FREQSUB = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load simulated data %%%

%%% W
CONDITION = 1;
SETSIGMA = 0.12;
c1 = 'W';
file = sprintf('metrics_PAR_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d_trapz.mat',...
                 CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
load(file)

MX_par_1_2(:,1) = MX_par';
Md_par_1_2(:,1) = Md_par';
Mi_par_1_2(:,1) = Mi_par';
xMX_par_1_2(:,1) = xMX_par';
xMd_par_1_2(:,1) = xMd_par';
xMi_par_1_2(:,1) = xMi_par';

%%% N3
CONDITION = 2;
SETSIGMA = 0.06;
c2 = 'N3';
file = sprintf('metrics_PAR_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d_trapz.mat',...
                 CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
load(file)

MX_par_1_2(:,2) = MX_par';
Md_par_1_2(:,2) = Md_par';
Mi_par_1_2(:,2) = Mi_par';
xMX_par_1_2(:,2) = xMX_par';
xMd_par_1_2(:,2) = xMd_par';
xMi_par_1_2(:,2) = xMi_par';

cond1 = cellstr(repmat(c1,90,1));
cond2 = cellstr(repmat(c2,90,1));

dataMX = cat(1,MX_par_1_2(:,1),MX_par_1_2(:,2));
dataMd = cat(1,Md_par_1_2(:,1),Md_par_1_2(:,2));
dataMi = cat(1,Mi_par_1_2(:,1),Mi_par_1_2(:,2));
dataxMX = cat(1,xMX_par_1_2(:,1),xMX_par_1_2(:,2));
dataxMd = cat(1,xMd_par_1_2(:,1),xMd_par_1_2(:,2));
dataxMi = cat(1,xMi_par_1_2(:,1),xMi_par_1_2(:,2));
condis = cat(1,cond1,cond2);


%%%% FIGURES %%%%
%
clf
%
%%% Integral xMi
% VIOLIN
violinplot(dataxMi,condis,'ShowMean',true), xlabel('Condition'), ylabel('Integral Violation of FDT')
legend('','','',num2str(mean(xMi_par_1_2(:,2)),'%.3e'),'','','','','','','',num2str(mean(xMi_par_1_2(:,1)),'%.3e'),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMi_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMi_viol.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(xMi_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Integral Violation of FDT'), hold on
plot(xMi_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(xMi_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(xMi_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(xMi_par_1_2(:,2)),'%.3e'),num2str(mean(xMi_par_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMi_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMi_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%% Differential xMd
% VIOLIN
violinplot(dataxMd,condis,'ShowMean',true), xlabel('Condition'), ylabel('Differential Violation of FDT')
legend('','','',num2str(mean(xMd_par_1_2(:,2)),'%.3e'),'','','','','','','',num2str(mean(xMd_par_1_2(:,1)),'%.3e'),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMd_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMd_viol.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(xMd_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Differential Violation of FDT'), hold on
plot(xMd_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(xMd_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(xMd_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(xMd_par_1_2(:,2)),'%.3e'),num2str(mean(xMd_par_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMd_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMd_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%% Integral Mi
% VIOLIN
violinplot(dataMi,condis,'ShowMean',true), xlabel('Condition'), ylabel('Integral Violation of FDT')
legend('','','',num2str(mean(Mi_par_1_2(:,2)),'%.3e'),'','','','','','','',num2str(mean(Mi_par_1_2(:,1)),'%.3e'),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Mi_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Mi_viol.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(Mi_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Integral Violation of FDT'), hold on
plot(Mi_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(Mi_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(Mi_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(Mi_par_1_2(:,2)),'%.3e'),num2str(mean(Mi_par_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Mi_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Mi_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%% Differential Md
% VIOLIN
violinplot(dataMd,condis,'ShowMean',true), xlabel('Condition'), ylabel('Differential Violation of FDT')
legend('','','',num2str(mean(Md_par_1_2(:,2)),'%.3e'),'','','','','','','',num2str(mean(Md_par_1_2(:,1)),'%.3e'),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Md_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Md_viol.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(Md_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Differential Violation of FDT'), hold on
plot(Md_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(Md_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(Md_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(Md_par_1_2(:,2)),'%.3e'),num2str(mean(Md_par_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Md_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Md_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%% FDR
% VIOLIN
violinplot(dataMX,condis,'ShowMean',true), xlabel('Condition'), ylabel('Fluctuation Dissipation Ratio')
legend('','','',num2str(mean(MX_par_1_2(:,2))),'','','','','','','',num2str(mean(MX_par_1_2(:,1))),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_MX_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_MX_viol.png',NSUBSIM);
saveas(f,file)
close(f)
% PLOT
plot(MX_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Differential Violation of FDT'), hold on
plot(MX_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(MX_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(MX_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(MX_par_1_2(:,2))),num2str(mean(MX_par_1_2(:,1))),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_MX_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_MX_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)
%%% xFDR
% VIOLIN
violinplot(dataxMX,condis,'ShowMean',true), xlabel('Condition'), ylabel('Fluctuation Dissipation Ratio')
legend('','','',num2str(mean(xMX_par_1_2(:,2))),'','','','','','','',num2str(mean(xMX_par_1_2(:,1))),'','','','Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMX_viol.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMX_viol.png',NSUBSIM);
saveas(f,file)
close(f)
% PLOT
plot(xMX_par_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Differential Violation of FDT'), hold on
plot(xMX_par_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(xMX_par_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(xMX_par_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(xMX_par_1_2(:,2))),num2str(mean(xMX_par_1_2(:,1))),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_xMX_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_xMX_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%%%%%%%%%%%%%%% CREATE DATAFILES FOR DOING VIOLIN PLOTS IN PYTHON
%%% xMi
file = sprintf('Sleep_NSIM_%d_xMi_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},xMi_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},xMi_par_1_2(i,2)); end;
fclose(fileID);
%%% xMd
file = sprintf('Sleep_NSIM_%d_xMd_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},xMd_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},xMd_par_1_2(i,2)); end;
fclose(fileID);
%%% xMX
file = sprintf('Sleep_NSIM_%d_xMX_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},xMX_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},xMX_par_1_2(i,2)); end;
fclose(fileID);
%%% Mi
file = sprintf('Sleep_NSIM_%d_Mi_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},Mi_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},Mi_par_1_2(i,2)); end;
fclose(fileID);
%%% Md
file = sprintf('Sleep_NSIM_%d_Md_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},Md_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},Md_par_1_2(i,2)); end;
fclose(fileID);
%%% MX
file = sprintf('Sleep_NSIM_%d_MX_par_for_violin.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:15, fprintf(fileID,'%s,%f\n',cond1{i},MX_par_1_2(i,1)); end;
for i = 1:15, fprintf(fileID,'%s,%f\n',cond2{i},MX_par_1_2(i,2)); end;
fclose(fileID);

end %NSUBSIM


%%% STATISTICAL TEST

% [h, p] = ttest2(xMd_par_1_2(:,1),xMd_par_1_2(:,2),'Alpha',0.05,'Vartype','unequal')
% [p, h] = ranksum(xMi_par_1_2(:,1),xMi_par_1_2(:,2))

%%%%%%%%%%%%%%%%%%%%%%
