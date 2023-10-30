clear all;
%% Load elementary effects data
date2save = "2023-10-30"; %"2023-10-27";
notes = "SS1"; %'MASS1';

% plas_conc
var = "plas_conc"
fname = strcat("./results_final/", date2save, '_MorrisSS_EE',...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tplas = readtable(fname, 'ReadRowNames', true);
EEplas = table2array(Tplas);

% musc_conc
var = "mus_conc"
fname = strcat("./results_final/", date2save, '_MorrisSS_EE',...
            '_var-', var,...
            "_notes-", notes, ".csv");
Tmusc = readtable(fname, 'ReadRowNames', true);
EEmusc = table2array(Tmusc);

parnames = Tmusc.Properties.VariableNames;
pname_plt = {};
for ii = 1:length(parnames)
    pname_plt{ii} = change_parname(parnames{ii});
end
%% Preprocess elementary effects data
% remove elementary effects above a threshold for computation problems
th = 900;
[r,c] = find(abs(EEplas) > th);
rplas_unique = unique(r);
EEplas(rplas_unique, :) = 0;
fprintf('removed %i values from EEplas \n', length(rplas_unique))

[r,c] = find(abs(EEmusc) > th);
rmusc_unique = unique(r);
EEmusc(rmusc_unique, :) = 0;
fprintf('removed %i values from EEmusc \n', length(rmusc_unique))

%% compute mu, mustar, sigma for each parameter
% row 1: mu, row 2: mustar, row3: sigma, row4: MI
plasvals = zeros(4, length(parnames));
muscvals = zeros(4, length(parnames));

plasvals(1,:) = mean(EEplas);
plasvals(2,:) = mean(abs(EEplas));
plasvals(3,:) = std(EEplas);
plasvals(4,:) = sqrt(plasvals(1,:).^2 + plasvals(3,:).^2);

muscvals(1,:) = mean(EEmusc);
muscvals(2,:) = mean(abs(EEmusc));
muscvals(3,:) = std(EEmusc);
muscvals(4,:) = sqrt(muscvals(1,:).^2 + muscvals(3,:).^2);

% IDs for significant MI
th = 0.1; % threshold for "significance"
temp = find(plasvals(4,:) > 0.1);
temp2 = find(muscvals(4,:) > 0.1);

MIsig = union(temp, temp2);

%% make figures
% Morrsi Plots
cmap = turbo(length(MIsig));
marksize=10; ms = 'diamond';
fx = 35; fy = 35; fleg = 16; ft = 22;
ftxt = 16; fgca = 18;
figure(1)
clf
nr = 1; nc = 2;
% plas_con
subplot(nr, nc, 1)
hold on
for jj = 1:length(MIsig)
    id = MIsig(jj);
    plot(plasvals(2, id), plasvals(3, id),...
            'markersize', marksize, 'marker', ms,...
            'color', cmap(jj,:), 'MarkerFaceColor', cmap(jj,:),...
            'linestyle', 'none', 'DisplayName', pname_plt{id})
    if or(plasvals(2,id) > 0.5, plasvals(3,id) > 0.1)
        if id == 12 % AdtKsec
            dx = -0.2; dy = 0.12;
        elseif id == 22 % ALD_eq
            dx = -0.25; dy = 0.1;
        elseif id == 15 % AcdKsec
            dx = -0.25; dy = 0.12;
        elseif id == 17 % AcdKreab
            dx = -0.15; dy = -0.1;
        elseif id == 16 % BcdKsec
            dx = 0.05; dy = 0.07;
        elseif id == 7 % mKaldo
            dx = -0.15; dy = -0.15;
        elseif id == 21 % Kecf_base
            dx = -0.1; dy = -0.12;
        elseif id == 11 % phi_dtKsec_eq
            dx = -0.1; dy = 0.15;
        elseif id == 10 % phi_GFR
            dx = -0.1; dy = 0.12;
        elseif plasvals(2,id) > 2.0
            dx = -0.22; dy = -0.045;
        else
            dx = -0.1; dy = 0.125;
        end
        text(plasvals(2,id) + dx, plasvals(3,id) + dy, pname_plt{id},...
            'fontsize', ftxt)
    end
end
set(gca, 'fontsize', fgca)
xlabel('\mu*', 'fontsize', fx)
ylabel('\sigma', 'fontsize', fy)
xlim([0.0, 3.0])
ylim([0.0, 3.0])
title('Plasma [K^+]', 'fontsize', ft)
grid on

% mus_con
subplot(nr, nc, 2)
hold on
for jj = 1:length(MIsig)
    id = MIsig(jj);
    plot(muscvals(2, id), muscvals(3, id),...
            'markersize', marksize, 'marker', ms,...
            'color', cmap(jj,:), 'MarkerFaceColor', cmap(jj,:),...
            'linestyle', 'none', 'DisplayName', pname_plt{id})
    if or(muscvals(2,id) > 20, muscvals(3,id) > 20)
        if id == 12 % AdtKsec
            dx = 4; dy = 0;
        elseif id == 22 % ALD_eq
            dx = -2; dy = 7;
        elseif id == 17 % AcdKreab
            dx = -3; dy = 5;
        elseif id == 16 % BcdKsec
            dx = 0; dy = -7;
        elseif id == 15 % AcdKsec
            dx = -10; dy = -3;
        elseif muscvals(2,id) > 60
            dx = -2; dy = -7;
        elseif muscvals(3,id) > 20
            dx = 4; dy = 0;
        else
            dx = 4; dy = 0;
        end
        text(muscvals(2,id) + dx, muscvals(3,id) + dy, pname_plt{id},...
            'fontsize', ftxt)
    end
end
set(gca, 'fontsize', fgca)
xlabel('\mu*', 'fontsize', fx)
ylabel('\sigma', 'fontsize', fy)
ylim([0.0, 150])
xlim([0.0, 150])
xticks(0.0:25:150)
yticks(0.0:25:150)
title('Intracellular [K^+]', 'fontsize', ft)
grid on

legend('fontsize',fleg)
AddLetters2Plots(figure(1), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)

%% Morris Index Plots
figure(2)
clf
nr = 2; nc = 1;
cmap = turbo(18);
marksize=15; ms = 'o';
fx = 22; fy = 22; fleg = 16; ft = 22;
ftxt = 16; fgca = 18;
cid = 16; cid2 = cid;
% plas_conc
subplot(nr,nc,1)
[~, Mids] = sort(plasvals(4,:), 'descend');
for ii = length(Mids):-1:1
    if ~ismember(Mids(ii), MIsig)
        Mids(ii) = [];
    end
end
plot(plasvals(4,Mids), 'linestyle', 'none', 'markersize', marksize,...
            'marker', ms, 'color', cmap(cid,:),...
            'markerfacecolor', cmap(cid,:))
set(gca, 'fontsize', 16)
xticks(1:length(Mids))
xlim([1,length(Mids)])
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
title('Plasma [K^+]', 'fontsize', ft)
xticklabels(pname_plt(Mids))
ylim([0.0,3.5])
grid on


% musc_conc
subplot(nr,nc,2)
plot(muscvals(4,Mids), 'linestyle', 'none','markersize',marksize,...
            'marker', ms, 'color', cmap(cid2,:),...
            'markerfacecolor', cmap(cid2,:))
set(gca, 'fontsize', 16)
xticks(1:length(Mids))
xlim([1,length(Mids)])
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
title('Intracellular [K^+]', 'fontsize', ft)
xticklabels(pname_plt(Mids))
grid on


AddLetters2Plots(figure(2), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)