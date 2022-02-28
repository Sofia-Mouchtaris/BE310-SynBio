% import true trial data
data = readtable('synbio_data.xlsx');

global conc

% normalize data to max and calculate standard error of mean
conc = data{:,1} * 20/2020;
gfpNorm = data{:,2:end} ./ max(data{:,2:end});
gfpSEM1 = std([gfpNorm(:,1), gfpNorm(:,4)],0,2) / sqrt(2);
gfpSEM2 = std([gfpNorm(:,2), gfpNorm(:,5)],0,2) / sqrt(2);
gfpSEM3 = std([gfpNorm(:,3), gfpNorm(:,6)],0,2) / sqrt(2);

%% define vars

global ror DELr n1 Atxgfp DELtxgfp Agfp DELgfp
ror = 0.5; %  μM^-3 min^-1
% luxr = .1; % uM
DELr = .0231; % δR = 0.0231 min-1
% Kr = 1.3 * 10^-5; % KR = 1.3e-5 μM
n1 = 1;
Atxgfp = .05; % αTXGFP = 0.05 μM/min
DELtxgfp = .2; % δTXGFP = 0.2 min-1
Agfp = 2; % αGFP = 2 min-1
DELgfp = 4 * 10^-4; % δGFP = 4e-4 min-1

%% optimize parameters

gfp = (data{:,2:4} + data{:,5:end}) ./ 2;
gfp = gfp ./ max(gfp);

currMin = [0,0,0,0];
currCost = inf;

% define range for each of 4 variables; loop over all possible combinations

for lowLux = linspace(.0004, .0006, 10)
for highLux = linspace(.13, .2, 10)
for lowKR = linspace(.05*10^-6, .3*10^-6, 10)
for highKR = linspace(1*10^-5, 3*10^-5, 10)

% calculate model

s1pr = [highLux, highKR];
s2pr = [lowLux, highKR];
s3pr = [lowLux, lowKR];

s1mod = Model(s1pr);
s2mod = Model(s2pr);
s3mod = Model(s3pr);

% calculate RMSE

s1cost = RMSE(s1mod, gfp(:,1));
s2cost = RMSE(s2mod, gfp(:,2));
s3cost = RMSE(s3mod, gfp(:,3));

totalCost = sqrt((s1cost^2 + s2cost^2 + s3cost^2)/3);

if totalCost < currCost
    currMin = [lowLux, highLux, lowKR, highKR];
    currCost = totalCost;
end

end
end
end
end

%% plot final model

conc = data{:,1} * 20/2020;

lowLux = 4.67 * 10^-4;
highLux = 1.69 * 10^-1;
lowKR  = 1.06 * 10^-7;
highKR = 2.11 * 10^-5;

s1pr = [highLux, highKR];
s2pr = [lowLux, highKR];
s3pr = [lowLux, lowKR];

s1mod = Model(s1pr);
s2mod = Model(s2pr);
s3mod = Model(s3pr);

s1cost = RMSE(s1mod, gfp(:,1));
s2cost = RMSE(s2mod, gfp(:,2));
s3cost = RMSE(s3mod, gfp(:,3));

oldconc = data{:,1} * 20/2020;
conc = logspace(-4, 4, 100) * 20/2020;

s1mod = Model(s1pr);
s2mod = Model(s2pr);
s3mod = Model(s3pr);

figure(1)
tiledlayout(3,1)
nexttile
semilogx(oldconc, gfp(:,1), 'b', 'LineWidth', 2)
hold on
errorbar(oldconc, gfp(:,1), gfpSEM1, 'b', 'LineWidth', 1.5)
semilogx(conc, s1mod, 'r', 'LineWidth', 2)
hold off
xlabel('[AHL] (uM)', 'FontSize', 22)
ylabel('Normalized [GFP]', 'FontSize', 22)
legend({'Data','', 'Model'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 1 - RMSE : %.3f', s1cost), 'FontSize', 24)

nexttile
semilogx(oldconc, gfp(:,2), 'b', 'LineWidth', 2)
hold on
errorbar(oldconc, gfp(:,2), gfpSEM2, 'b', 'LineWidth', 1.5)
semilogx(conc, s2mod, 'r', 'LineWidth', 2)
hold off
xlabel('[AHL] (uM)', 'FontSize', 22)
ylabel('Normalized [GFP]', 'FontSize', 22)
legend({'Data','','Model'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 2 - RMSE : %.3f', s2cost), 'FontSize', 24)

nexttile
semilogx(oldconc, gfp(:,3), 'b', 'LineWidth', 2)
hold on
errorbar(oldconc, gfp(:,3), gfpSEM3, 'b', 'LineWidth', 1.5)
semilogx(conc, s3mod, 'r', 'LineWidth', 2)
hold off
xlabel('[AHL] (uM)', 'FontSize', 22)
ylabel('Normalized [GFP]', 'FontSize', 22)
legend({'Data','', 'Model'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 3 - RMSE : %.3f', s3cost), 'FontSize', 24)


%% FUNCTIONS

function Yexp = Model(feat)
global conc ror DELr n1 Atxgfp DELtxgfp Agfp DELgfp
Yexp = zeros(length(conc),1); % calculate expected model values

for i_c = 1:length(conc)
Yexp(i_c) = Agfp / (DELtxgfp * DELgfp) * (Atxgfp * (ror * feat(1)^2 * conc(i_c)^2 / DELr / feat(2))^n1)...
    / (1 + (ror * feat(1)^2 * conc(i_c)^2 / DELr / feat(2))^n1);
end

Yexp = Yexp ./ max(Yexp); % normalize to max val

end

function cost = RMSE(Yexp, Y)

cost = sqrt(sum((Yexp - Y).^2)/length(Y));

end