% KCl Only Morris Plots figures for write up
clear all;

%% Load data 
date2save = "2023-10-27";
notes = "MA1";
sim_type = "KClOnly";

% amt_gut
var = "amt_gut"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis', ...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tgut = readtable(fname, 'ReadRowNames',true);

% conc_plas
var = "conc_plas"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tplas = readtable(fname,'ReadRowNames',true);

% conc_inter
var = "conc_inter"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tint = readtable(fname,'ReadRowNames',true);

% conc_muscle
var = "conc_muscle"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tmusc = readtable(fname,'ReadRowNames',true);

% parameter names list
parnames = ["V_plasma", ...
            "V_inter", ...
            "V_muscle", ...
            "kgut", ...
            "Km", ... 
            "Vmax", ...
            "m_K_ALDO", ...
            "P_ECF", ...
            "FF", ...
            "GFR", ...
            "dtKsec_eq", ...
            "A_dtKsec", ...
            "B_dtKsec", ...
            "cdKsec_eq", ...
            "A_cdKsec", ...
            "B_cdKsec", ...
            "A_cdKreab", ...
            "A_insulin", ...
            "B_insulin", ...
            "KMuscleBase",...
            "Kecf_base",...
            "ALD_eq",...
            "etapsKreab"
            ];


%% get values
tvals_plas = get_time_vals(parnames, Tplas);
tvals_musc = get_time_vals(parnames, Tmusc);
tvals_inter = get_time_vals(parnames, Tint);
tvals_gut = get_time_vals(parnames, Tgut);

parnames_plt = cell(size(parnames));
for ii = 1:length(parnames)
    parnames_plt{ii} = change_parname(parnames(ii));
end

%% Time series plot, remove low parameters
times = Tplas("time", :);
tvals = (table2array(times) - 460)./60;
nrows = 2; ncols = 2;
% find parameters that have a minimum sigma and mustar value
plas_sigvals_all = tvals_plas(:,:,3);
plas_mustarvals_all = tvals_plas(:,:,2);
plas_muvals_all = tvals_plas(:,:,1);
musc_sigvals_all = tvals_musc(:,:,3);
musc_mustarvals_all = tvals_musc(:,:,2);
musc_muvals_all = tvals_musc(:,:,1);

sig_min = 0.5; mus_min = 0.5;

allIDs = [];
temp = max(plas_sigvals_all, [], 2);
psig_ids = find(temp > sig_min);
allIDs = union(allIDs, psig_ids);
temp = max(plas_mustarvals_all, [], 2);
pmus_ids = find(temp > mus_min);
allIDs = union(allIDs, pmus_ids);

temp = max(musc_sigvals_all, [], 2);
msig_ids = find(temp > sig_min);
allIDs = union(allIDs, msig_ids);
temp = max(musc_mustarvals_all, [], 2);
mmus_ids = find(temp > mus_min);
allIDs = union(allIDs, mmus_ids);

% make plot
figure(1)
clf
cmap = turbo(length(allIDs));
xminmax = [-6, 8];
yminmax = [0.0, 3.2];
 ms = '.';
fx = 22; fy = 35; fleg = 16; ft = 22;
ftxt = 16; fgca = 18;
lw = 3; ls = {'-', '--', ':'}; % linestyles
% K plas mu star
subplot(nrows,ncols,1)
hold on
plt_parnames_sub = [];
for jj = 1:length(allIDs)
    ii = allIDs(jj);
    mustarvals = tvals_plas(ii, :, 2);
    lstyle = mod(jj,3) + 1;
    plot(tvals, mustarvals, 'linewidth', lw,...
        'color', cmap(jj,:),...
        'linestyle', ls{lstyle})
    plt_parnames_sub = [plt_parnames_sub; change_parname(parnames(ii))];
end
set(gca, 'fontsize', fgca)
xlabel('time', 'fontsize',fx)
ylabel('\mu*', 'fontsize',fy)
title("Plasma [K^+]", 'fontsize',ft)
xlim(xminmax)
ylim(yminmax)
%yticks(ytickvals)
grid on


% K muscle mu star
subplot(nrows, ncols,2)
hold on
for jj = 1:length(allIDs)
    ii = allIDs(jj);
    mustarvals = tvals_musc(ii, :, 2);
    lstyle = mod(jj,3) + 1;
    plot(tvals, mustarvals, 'linewidth', lw,...
        'color', cmap(jj,:),...
        'linestyle', ls{lstyle})
end
set(gca, 'fontsize', fgca)
xlabel('time', 'fontsize',fx)
xlim(xminmax)
ylim(yminmax)
%yticks(ytickvals)
ylabel('\mu*', 'fontsize',fy)
title("Intracellular [K^+]", 'fontsize',ft)
grid on


% time versus sigma


% K plas sigma
subplot(nrows, ncols,3)
hold on
for jj = 1:length(allIDs)
    ii = allIDs(jj);
    sigvals = tvals_plas(ii, :, 3);
    lstyle = mod(jj,3) + 1;
    plot(tvals, sigvals, 'linewidth', lw,...
        'color', cmap(jj,:), 'linestyle', ls{lstyle})
end
set(gca, 'fontsize', fgca)
xlabel('time', 'fontsize', fx)
xlim(xminmax)
ylim(yminmax)
%yticks(ytickvals)
ylabel('\sigma', 'fontsize',fy)
title("Plasma [K^+]", 'fontsize',ft)
grid on


% K muscle sigma
subplot(nrows, ncols,4)
hold on
for jj = 1:length(allIDs)
    ii = allIDs(jj);
    sigvals = tvals_musc(ii, :, 3);
    lstyle = mod(jj,3) + 1;
    plot(tvals, sigvals, 'linewidth', lw,...
        'color', cmap(jj,:), 'linestyle', ls{lstyle})
end
set(gca, 'fontsize', fgca)
xlabel('time', 'fontsize', fx)
xlim(xminmax)
ylim(yminmax)
%yticks(ytickvals)
ylabel('\sigma', 'fontsize',fy)
title("Intracellular [K^+]", 'fontsize',ft)
grid on

legend(plt_parnames_sub, 'fontsize', fleg)
%sgtitle("Meal + KCl")

AddLetters2Plots(figure(1), {'(A1)', '(B1)', '(A2)', '(B2)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)

%% Morris Indices
npars = length(parnames);
ntimes = size(times, 2);
plas_MIvals_all = zeros(npars, ntimes);
musc_MIvals_all = zeros(npars, ntimes);

for ii = 1:npars
    for jj = 1:ntimes
        plas_MIvals_all(ii, jj) = MorrisIndex(plas_muvals_all(ii, jj), plas_sigvals_all(ii,jj));
        musc_MIvals_all(ii, jj) = MorrisIndex(musc_muvals_all(ii,jj), musc_sigvals_all(ii,jj));
    end
end

plas_maxMI = max(plas_MIvals_all, [], 2);
[pMI, pIDs] = sort(plas_maxMI, "descend"); % order MI values
pars_plas_sort = [];
for ii = 1:length(pIDs)
    pID = pIDs(ii);
    pars_plas_sort = [pars_plas_sort; change_parname(parnames(pID))];
end

plas_meanMI = mean(plas_MIvals_all, 2);
plas_sdMI = std(plas_MIvals_all, 0, 2);
plas_minMI = min(plas_MIvals_all, [], 2);



musc_maxMI = max(musc_MIvals_all, [], 2);
[mMI, mIDs] = sort(musc_maxMI, "descend"); % order MI values
pars_musc_sort = [];
for ii = 1:length(mIDs)
    mID = mIDs(ii);
    pars_musc_sort = [pars_musc_sort; change_parname(parnames(mID))];
end
musc_meanMI = mean(musc_MIvals_all, 2);
musc_sdMI = std(musc_MIvals_all, 0, 2);
musc_minMI = min(musc_MIvals_all, [], 2);



% Morris Index plot
figure(2)
marksize = 15;
marksize2 = 12;
cmap = turbo(18);
cmax = cmap(16,:);%cmap(3,:);
cmean = cmap(16,:);
yminmax = [0.0, 4.25];
ms1 = '*'; % for minmax
ms2 = 'o'; % for mean
nr = 2; nc = 1;
clf;
subplot(nr, nc, 1)
hold on


plot(plas_meanMI(pIDs), 'markersize', marksize, 'marker', ms2,...
            'color', cmean, 'markerfacecolor', cmean,'linestyle', 'none')
errorbar(plas_meanMI(pIDs), plas_sdMI(pIDs),...
            'color', cmean, 'linestyle', 'none', 'linewidth',2.0)
plot(plas_maxMI(pIDs), 'markersize', marksize2, 'marker', ms1,...
            'color', cmax, 'markerfacecolor', cmax, 'linestyle', 'none')
plot(plas_minMI(pIDs), 'markersize', marksize2, 'marker', ms1,...
            'color', cmax, 'markerfacecolor',cmax,'linestyle', 'none')
xticks(1:23)
xticklabels(pars_plas_sort)
xlim([1,23])
title("Plasma [K^+]")
ylabel("Morris Index")
legend('mean', '', 'maximum/minimum')
set(gca, 'fontsize', 18)
ylim(yminmax)
%yticks(ytickvals)
grid on

subplot(nr, nc, 2)
hold on
plot(musc_meanMI(pIDs), 'markersize', marksize, 'marker', ms2,...
        'color', cmean, 'markerfacecolor',cmean, 'linestyle', 'none')
errorbar(musc_meanMI(pIDs), musc_sdMI(pIDs),...
        'color', cmean, 'linestyle', 'none','linewidth',2.0)
plot(musc_maxMI(pIDs), 'markersize', marksize2, 'marker', ms1,...
            'color', cmax,'markerfacecolor',cmax, 'linestyle', 'none')
plot(musc_minMI(pIDs), 'markersize', marksize2, 'marker', ms1,...
        'color', cmax, 'markerfacecolor',cmax,'linestyle', 'none')
xticks(1:23)
xlim([1,23])
xticklabels(pars_plas_sort)
title("Intracellular [K^+]")
legend('mean', '', 'maximum/minimum')
ylabel("Morris Index")
set(gca, 'fontsize', 18)
ylim(yminmax)
%yticks(ytickvals)
grid on


AddLetters2Plots(figure(2), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)


%----------------------
% functions used
%----------------------
%% Morris Index
function MI = MorrisIndex(mustar, sigma)
    MI = sqrt(mustar^2 + sigma^2);
end

%% get time vals
function tvals = get_time_vals(parnames, T)
    ntimes = size(T.Properties.VariableNames, 2);
    tvals = zeros(length(parnames), ntimes, 3); % 1: mu, 2: mu*, 3: sigma
    for ii = 1:length(parnames)
        pname = parnames(ii);
        tvals(ii,:,1) = get_mu(T, pname);
        tvals(ii,:,2) = get_mustar(T, pname);
        tvals(ii,:,3) = get_sigma(T, pname);
    end
end

%% get sigma 
function sig_vals = get_sigma(T, pname)
    % Input
    %   T - table
    %   pname - parameter name
    nm = strcat('sigma_', pname);
    nt = size(T.Properties.VariableNames, 2);
    vals = T(nm, :);
    sig_vals = zeros(nt,1);
    for ii = 1:nt
        temp = strcat('vals.time', num2str(ii));
        sig_vals(ii) = eval(temp);
    end
end

%% get mu star values
function mus_vals = get_mustar(T, pname)
    % Input
    %   T - table
    %   pname - parameter name
    nm = strcat('mu.star_', pname);
    nt = size(T.Properties.VariableNames, 2);
    vals = T(nm, :);
    mus_vals = zeros(nt,1);
    for ii = 1:nt
        temp = strcat('vals.time', num2str(ii));
        mus_vals(ii) = eval(temp);
    end
end

%% get mu values
function mu_vals = get_mu(T, pname)
    % Input
    %   T - table
    %   pname - parameter name
    nm = strcat('mu_', pname);
    nt = size(T.Properties.VariableNames, 2);
    vals = T(nm, :);
    mu_vals = zeros(nt, 1);
    for ii = 1:nt
        temp = strcat('vals.time', num2str(ii));
        mu_vals(ii) = eval(temp);
    end
end

