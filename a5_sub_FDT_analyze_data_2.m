%%% A5_SUB_FDT_ANALYZE_DATA_2  Plots subject-wise FDT metrics: Wakefulness vs Deep Sleep.
% Author: Juan Manuel Monti
%
% Loads metrics_SUB_*.mat for both conditions and generates boxplots and line
% plots comparing W vs N3 across subjects for Mi (integral), Md (differential),
% and MX (FDT ratio) metrics.
%
% Figures saved per metric (Mi, Md, MX):
%   Sleep_NSIM_<N>_<metric>_box.{svg,png}   boxplot (W vs N3)
%   Sleep_NSIM_<N>_<metric>_plot.{svg,png}  per-subject line plot
%
% CSV export (remove the `return` on line ~163 to activate):
%   Sleep_NSIM_<N>_<metric>_sub.csv  columns: cond, value
%
% Requires: metrics_SUB_Sleep_COND_*.mat

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
file = sprintf('metrics_SUB_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d.mat',...
                 CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
load(file)

MX_sub_1_2(:,1) = MX_sub';
Md_sub_1_2(:,1) = Md_sub';
Mi_sub_1_2(:,1) = Mi_sub';

xMX_sub_1_2(:,1) = xMX_sub';
xMd_sub_1_2(:,1) = xMd_sub';
xMi_sub_1_2(:,1) = xMi_sub';

%%% N3
CONDITION = 2;
SETSIGMA = 0.06;
c2 = 'N3';
file = sprintf('metrics_SUB_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d.mat',...
                 CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
load(file)

MX_sub_1_2(:,2) = MX_sub';
Md_sub_1_2(:,2) = Md_sub';
Mi_sub_1_2(:,2) = Mi_sub';

xMX_sub_1_2(:,2) = xMX_sub';
xMd_sub_1_2(:,2) = xMd_sub';
xMi_sub_1_2(:,2) = xMi_sub';

cond1 = cellstr(repmat(c1,NSUB,1));
cond2 = cellstr(repmat(c2,NSUB,1));

dataMX = cat(1,MX_sub_1_2(:,1),MX_sub_1_2(:,2));
dataMd = cat(1,Md_sub_1_2(:,1),Md_sub_1_2(:,2));
dataMi = cat(1,Mi_sub_1_2(:,1),Mi_sub_1_2(:,2));

dataxMX = cat(1,xMX_sub_1_2(:,1),xMX_sub_1_2(:,2));
dataxMd = cat(1,xMd_sub_1_2(:,1),xMd_sub_1_2(:,2));
dataxMi = cat(1,xMi_sub_1_2(:,1),xMi_sub_1_2(:,2));
condis = cat(1,cond1,cond2);

%%%% FIGURES %%%%
clf

%%% Integral Mi
% BOXPLOT
boxplot(dataMi,condis), xlabel('Condition'), ylabel('Integral Violation of FDT')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Mi_box.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Mi_box.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(Mi_sub_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Integral Violation of FDT'), hold on
plot(Mi_sub_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(Mi_sub_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(Mi_sub_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(Mi_sub_1_2(:,2)),'%.3e'),num2str(mean(Mi_sub_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Mi_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Mi_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)

%%% Differential Md
% BOXPLOT
boxplot(dataMd,condis), xlabel('Condition'), ylabel('Integral Violation of FDT')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Md_box.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Md_box.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(Md_sub_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Integral Violation of FDT'), hold on
plot(Md_sub_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(Md_sub_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(Md_sub_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(Md_sub_1_2(:,2)),'%.3e'),num2str(mean(Md_sub_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_Md_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_Md_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)


%%% Ratio MX
% BOXPLOT
boxplot(dataMX,condis), xlabel('Condition'), ylabel('Integral Violation of FDT')
f = gcf;
file = sprintf('Sleep_NSIM_%d_MX_box.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_MX_box.png',NSUBSIM);
saveas(f,file)
clf
close(f)
% PLOT
plot(MX_sub_1_2(:,2),'-s','LineWidth',1,'MarkerFaceColor',[.6 .6 1]), xlabel('Parcells'), ylabel('Integral Violation of FDT'), hold on
plot(MX_sub_1_2(:,1),'-s','LineWidth',1,'MarkerFaceColor',[1 .6  0])
yline(mean(MX_sub_1_2(:,2)),'LineWidth',1,'Color',[.6 .6 1])
yline(mean(MX_sub_1_2(:,1)),'LineWidth',1,'Color',[1 .6  0])
legend(num2str(mean(MX_sub_1_2(:,2)),'%.3e'),num2str(mean(MX_sub_1_2(:,1)),'%.3e'),'Location','NorthEast')
f = gcf;
file = sprintf('Sleep_NSIM_%d_MX_plot.svg',NSUBSIM);
print(f,'-dsvg',file)
file = sprintf('Sleep_NSIM_%d_MX_plot.png',NSUBSIM);
saveas(f,file)
clf
close(f)

return

%%% CREATE DATAFILES FOR DOING BOXPLOTS IN PYTHON %%%%%%%%%%%%%%%%%%%%
%%% xMi
file = sprintf('Sleep_NSIM_%d_xMi_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},xMi_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},xMi_sub_1_2(i,2)); end;
fclose(fileID);
%%% xMd
file = sprintf('Sleep_NSIM_%d_xMd_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},xMd_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},xMd_sub_1_2(i,2)); end;
fclose(fileID);
%%% MX
file = sprintf('Sleep_NSIM_%d_xMX_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},xMX_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},xMX_sub_1_2(i,2)); end;
fclose(fileID);
%%% Mi
file = sprintf('Sleep_NSIM_%d_Mi_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},Mi_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},Mi_sub_1_2(i,2)); end;
fclose(fileID);
%%% Md
file = sprintf('Sleep_NSIM_%d_Md_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},Md_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},Md_sub_1_2(i,2)); end;
fclose(fileID);
%%% MX
file = sprintf('Sleep_NSIM_%d_MX_sub.csv',NSUBSIM);
fileID = fopen(file,'w');
fprintf(fileID,'%s,%s\n','cond','value');
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond1{i},MX_sub_1_2(i,1)); end;
for i = 1:NSUB, fprintf(fileID,'%s,%f\n',cond2{i},MX_sub_1_2(i,2)); end;
fclose(fileID);

end % NSUBSIM
%%%%%%%%%%%%%%%%%%%%%%
return
