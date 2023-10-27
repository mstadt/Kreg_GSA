% Make plots for single meal simulations
clear all

%% Import simulation data
fprintf('loading data \n')
% Meal + KCl
fname1 = "2023-10-27_MealSim__type-MealKCl_notes-singlemeal.csv";
dat1   = readtable(strcat("./results_final/",fname1));
lab1   = "Meal + KCl";

% KCl only
fname2 = "2023-10-27_MealSim__type-KClOnly_notes-singlemeal.csv";
dat2   = readtable(strcat("./results_final/",fname2));
lab2   = 'KCl Only';

% Meal only
fname3 = "2023-10-27_MealSim__type-MealOnly_notes-singlemeal.csv";
dat3   = readtable(strcat("./results_final/",fname3));
lab3   = 'Meal Only';

%% Make figure
tshift = 460; % fasting time
% figure specs
lw = 3;
f.xlab = 18; f.ylab = 18; f.title = 20;
f.leg = 16; f.gca = 16;
cmap = parula(6);
c1 = cmap(1,:);
c2 = cmap(3,:);
c3 = cmap(4,:);
ls1 = '-'; ls2 = ':'; ls3 = '-.';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
x0 = -6; xf = 8;
leglabs = {lab1, lab2, lab3};

fprintf('making figs \n')


figure(1)
clf
nrows=2; ncols=2;
subplot(nrows,ncols,1)
hold on
plot((dat1.time - tshift)/60, dat1.amt_gut,'linewidth',lw,'color',c1,'linestyle',ls1)
plot((dat2.time - tshift)/60, dat2.amt_gut,'linewidth',lw,'color',c2,'linestyle',ls2)
plot((dat3.time - tshift)/60, dat3.amt_gut,'linewidth',lw,'color',c3,'linestyle',ls3)
xlim([x0,xf])
set(gca, 'fontsize', f.gca)
ylabel('M_{Kgut} (mmol)', 'fontsize',f.ylab)
xlabel('time (hours)','fontsize',f.xlab)
%title('Gut K^+ Amount', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

subplot(nrows,ncols,2)
hold on
plot((dat1.time - tshift)/60, dat1.conc_plas,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot((dat2.time - tshift)/60, dat2.conc_plas,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot((dat3.time - tshift)/60, dat3.conc_plas,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
ylabel('K_{plasma} (mmol/L)', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
%title('Plasma [K^+]', 'fontsize', f.title)
xlim([x0,xf])
grid on
legend(leglabs, 'fontsize', f.leg)

subplot(nrows,ncols,3)
hold on
plot((dat1.time - tshift)/60,dat1.conc_inter,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot((dat2.time - tshift)/60,dat2.conc_inter,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot((dat3.time - tshift)/60,dat3.conc_inter,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
ylabel('K_{inter} (mmol/L)', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
%title('Interstitial [K^+]', 'fontsize', f.title)
xlim([x0,xf])
grid on
legend(leglabs, 'fontsize', f.leg)

subplot(nrows,ncols,4)
hold on
plot((dat1.time - tshift)/60,dat1.conc_muscle,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot((dat2.time - tshift)/60,dat2.conc_muscle,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot((dat3.time - tshift)/60,dat3.conc_muscle,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
ylabel('K_{intracellular} (mmol/L)', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
%title('Intracellular [K^+]', 'fontsize', f.title)
xlim([x0,xf])
grid on
legend(leglabs, 'fontsize', f.leg)

AddLetters2Plots(figure(1), {'(A)', '(B)', '(C)', '(D)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 16)
