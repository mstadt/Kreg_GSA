% MI plot for all 3 meal types
clear all;

%% Load data
date2save = "2023-10-27";
notes = "MA1";

% Meal + KCl
sim_type = "MealKCl"

var = "conc_plas"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tplas_MealKCl = readtable(fname,'ReadRowNames',true);

var = "conc_muscle"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tmusc_MealKCl = readtable(fname,'ReadRowNames',true);

% Meal Only
sim_type = "MealOnly"

var = "conc_plas"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tplas_MealOnly = readtable(fname,'ReadRowNames',true);

var = "conc_muscle"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tmusc_MealOnly = readtable(fname,'ReadRowNames',true);

% KCl Only
sim_type = "KClOnly"

var = "conc_plas"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tplas_KClOnly = readtable(fname,'ReadRowNames',true);

var = "conc_muscle"
fname = strcat("./results_final/", ...
            date2save, '_MorrisAnalysis',...
            "_type-", sim_type, ...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tmusc_KClOnly = readtable(fname,'ReadRowNames',true);

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


%% Preprocess data
% Meal + KCl
[plas_MIvals, musc_MIvals] = getMIs(Tplas_MealKCl, Tmusc_MealKCl, parnames);
MIstats_plasMealKCl = getMI_stats(plas_MIvals);
MIstats_muscMealKCl = getMI_stats(musc_MIvals);

% Meal Only
[plas_MIvals, musc_MIvals] = getMIs(Tplas_MealOnly, Tmusc_MealOnly, parnames);
MIstats_plasMealOnly = getMI_stats(plas_MIvals);
MIstats_muscMealOnly = getMI_stats(musc_MIvals);

% KCl Only
[plas_MIvals, musc_MIvals] = getMIs(Tplas_KClOnly, Tmusc_KClOnly, parnames);
MIstats_plasKClOnly = getMI_stats(plas_MIvals);
MIstats_muscKClOnly = getMI_stats(musc_MIvals);

% Compute mean MI for all
meanMIs = [MIstats_plasMealKCl(:,1), MIstats_plasMealOnly(:,1), MIstats_plasKClOnly(:,1)];
meanMI_all = mean(meanMIs,2);
[sorted_MI, IDs] = sort(meanMI_all, 'descend');

parnames_sort = parnames(IDs);
parnames_plt_sort = cell(size(parnames_sort));
for ii = 1:length(parnames_sort)
    parnames_plt_sort{ii} = change_parname(parnames_sort(ii));
end

%% Make figure
figure(1)
clf
marksize = 15;
fleg = 18; fx = 18; fy = 18;
cmap = parula(4);
c1 = cmap(1,:);
c2 = cmap(2,:);
c3 = cmap(3,:);
ms1 = 'o';
ms2 = 'diamond';
ms3 = 'square';
nr = 2; nc = 1;

subplot(nr,nc,1)
hold on
plot(MIstats_plasMealKCl(IDs,1), 'markersize', marksize, 'marker', ms1,...
                                    'color', c1, 'markerfacecolor', c1, ...
                                    'linestyle', 'none')
errorbar(MIstats_plasMealKCl(IDs,1), MIstats_plasMealKCl(IDs,2), ...
                                    'color', c1, 'linestyle', 'none', 'linewidth', 2.0)
plot(MIstats_plasKClOnly(IDs,1), 'markersize', marksize, 'marker', ms2,...
                                    'color', c2, 'markerfacecolor', c2, ...
                                    'linestyle', 'none')
errorbar(MIstats_plasKClOnly(IDs,1), MIstats_plasKClOnly(IDs,2), ...
                                    'color', c2, 'linestyle', 'none', 'linewidth', 2.0)
plot(MIstats_plasMealOnly(IDs,1), 'markersize', marksize, 'marker', ms3,...
                                    'color', c3, 'markerfacecolor', c3, ...
                                    'linestyle', 'none')
errorbar(MIstats_plasMealOnly(IDs,1), MIstats_plasMealOnly(IDs,2), ...
                                    'color', c3, 'linestyle', 'none', 'linewidth', 2.0)
legend('Meal + KCl', '', ...
    'KCl Only', '', ...
    'Meal Only', '',...
    'fontsize', fleg)
set(gca, 'fontsize', 18)
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
ylim([0.0,4.5])
xticks(1:23)
xticklabels(parnames_plt_sort)
xlim([1,23])
grid on

subplot(nr,nc,2)
hold on
plot(MIstats_muscMealKCl(IDs,1), 'markersize', marksize, 'marker', ms1,...
                                    'color', c1, 'markerfacecolor', c1, ...
                                    'linestyle', 'none')
errorbar(MIstats_muscMealKCl(IDs,1), MIstats_muscMealKCl(IDs,2), ...
                                    'color', c2, 'linestyle', 'none', 'linewidth', 2.0)
plot(MIstats_muscKClOnly(IDs,1), 'markersize', marksize, 'marker', ms2,...
                                    'color', c2, 'markerfacecolor', c2, ...
                                    'linestyle', 'none')
errorbar(MIstats_muscKClOnly(IDs,1), MIstats_muscKClOnly(IDs,2), ...
                                    'color', c2, 'linestyle', 'none', 'linewidth', 2.0)
plot(MIstats_muscMealOnly(IDs,1), 'markersize', marksize, 'marker', ms3,...
                                    'color', c3, 'markerfacecolor', c3, ...
                                    'linestyle', 'none')
errorbar(MIstats_muscMealOnly(IDs,1), MIstats_muscMealOnly(IDs,2), ...
                                    'color', c3, 'linestyle', 'none', 'linewidth', 2.0)

legend('Meal + KCl', '', ...
    'KCl Only', '', ...
    'Meal Only', '',...
    'fontsize', fleg)
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
set(gca, 'fontsize', 18)
xticks(1:23)
xticklabels(parnames_plt_sort)
xlim([1,23])
ylim([0.0,4.5])
grid on

AddLetters2Plots(figure(1), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)
%%
% Functions used
function  MIstats = getMI_stats(MIvals)
    meanMI = mean(MIvals, 2);
    sdMI = std(MIvals, 0, 2);
    minMI = min(MIvals, [], 2);
    maxMI = max(MIvals, [], 2);
    MIstats = [meanMI, sdMI, minMI, maxMI];
end

function [plas_MIvals, musc_MIvals] = getMIs(Tplas, Tmusc, parnames)
    tvals_plas = get_time_vals(parnames, Tplas);
    tvals_musc = get_time_vals(parnames, Tmusc);
    plas_sigvals_all = tvals_plas(:,:,2);
    plas_muvals_all = tvals_plas(:,:,1);
    musc_sigvals_all = tvals_musc(:,:,2);
    musc_muvals_all = tvals_musc(:,:,1);
    
    times = Tplas("time", :);
    npars = length(parnames);
    ntimes = size(times, 2);
    
    plas_MIvals = zeros(npars, ntimes);
    musc_MIvals = zeros(npars, ntimes);
    for ii = 1:npars
        for jj = 1:ntimes
            plas_MIvals(ii, jj) = MorrisIndex(plas_muvals_all(ii, jj), plas_sigvals_all(ii,jj));
            musc_MIvals(ii, jj) = MorrisIndex(musc_muvals_all(ii,jj), musc_sigvals_all(ii,jj));
        end
    end
end

% Morris Index
function MI = MorrisIndex(mustar, sigma)
    MI = sqrt(mustar^2 + sigma^2);
end

% get time vals
function tvals = get_time_vals(parnames, T)
    ntimes = size(T.Properties.VariableNames, 2);
    tvals = zeros(length(parnames), ntimes, 2); % 1: mu, 2: sigma
    for ii = 1:length(parnames)
        pname = parnames(ii);
        tvals(ii,:,1) = get_mu(T, pname);
        tvals(ii,:,2) = get_sigma(T, pname);
    end
end

%get sigma 
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

% get mu values
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