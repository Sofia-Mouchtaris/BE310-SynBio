% define variables

global conc ror DELr n1 Atxgfp DELtxgfp Agfp DELgfp
ror = 0.5; %  μM^-3 min^-1
% luxr = .1; % uM
DELr = .0231; % δR = 0.0231 min-1
% Kr = 1.3 * 10^-5; % KR = 1.3e-5 μM
n1 = 1;
Atxgfp = .05; % αTXGFP = 0.05 μM/min
DELtxgfp = .2; % δTXGFP = 0.2 min-1
Agfp = 2; % αGFP = 2 min-1
DELgfp = 4 * 10^-4; % δGFP = 4e-4 min-1
conc = logspace(-4, 4, 100) * 20/2020;

lowLux = 4.67 * 10^-4;
highLux = 1.69 * 10^-1;
lowKR  = 1.06 * 10^-7;
highKR = 2.11 * 10^-5;

% create expected steady state models

s1pr = [highLux, highKR];
s2pr = [lowLux, highKR];
s3pr = [lowLux, lowKR];

s1mod = Model(s1pr);
s2mod = Model(s2pr);
s3mod = Model(s3pr);

[~, i] = max(diff(s1mod));
s1switch = conc(i);
[~, i] = max(diff(s2mod));
s2switch = conc(i);
[~, i] = max(diff(s3mod));
s3switch = conc(i);

%%

% define plate and model parameters

Radius_Plate = 40; % mm
Radius_Disk = 3; % mm
T = 24 * 60; % Experiment duration (minutes) 
D = .06; % Diffusion rate (mm^2/min)
sourceconc = 10; % Source concentration (uM) 10 μL of 10 μM AHL
dAHL = 4.8135 * 10^-6; % min-1 

meshpts = 101; 
dx = Radius_Plate / ((meshpts - 1)/2); 
dy = dx;
xgrid = -Radius_Plate:dx:Radius_Plate;
ygrid = -Radius_Plate:dy:Radius_Plate;
[X,Y] = meshgrid(xgrid, ygrid); 
mid = (meshpts-1)/2 + 1;

r = .25;
dt = r * (dx^2)/D;

AHL_Initial = zeros(length(xgrid), length(ygrid));
AHL_Initial(sqrt(X.^2+Y.^2)<=Radius_Disk) = sourceconc;

time = 0:dt:T;

% create empty x by y by t for all relevant variables

AHL = zeros(length(xgrid), length(ygrid), length(time));
R1 = zeros(length(xgrid), length(ygrid), length(time));
TX1 = zeros(length(xgrid), length(ygrid), length(time));
GFP1 = zeros(length(xgrid), length(ygrid), length(time));

R2 = zeros(length(xgrid), length(ygrid), length(time));
TX2 = zeros(length(xgrid), length(ygrid), length(time));
GFP2 = zeros(length(xgrid), length(ygrid), length(time));

R3 = zeros(length(xgrid), length(ygrid), length(time));
TX3 = zeros(length(xgrid), length(ygrid), length(time));
GFP3 = zeros(length(xgrid), length(ygrid), length(time));

EDGEAHL1 = zeros(1, length(time));
EDGEAHL2 =zeros(1, length(time));
EDGEAHL3 = zeros(1, length(time));

EDGEGFP1 = zeros(1, length(time));
EDGEGFP2 = zeros(1, length(time));
EDGEGFP3 = zeros(1, length(time));

% set initial conditions

AHL(:,:,1) = AHL_Initial;
EDGEAHL1(1) = getSwitch(AHL(mid,:,1), s1switch, xgrid);
EDGEAHL2(1) = getSwitch(AHL(mid,:,1), s2switch, xgrid);
EDGEAHL3(1) = getSwitch(AHL(mid,:,1), s3switch, xgrid);

gfpswitch = .1;

%%

for i_t = 2:length(time)
for i_x = 2:(length(xgrid)-1)
for i_y = 2:(length(ygrid)-1)

% calculate AHL diffusion

AHL(i_x, i_y, i_t) = AHL(i_x, i_y, i_t - 1) + r * (AHL(i_x - 1, i_y, i_t - 1) + ...
    AHL(i_x + 1, i_y, i_t - 1) + AHL(i_x, i_y - 1, i_t - 1) + AHL(i_x, i_y + 1, i_t - 1) - ...
    4 * AHL(i_x, i_y, i_t - 1)) - dt * (dAHL * AHL(i_x, i_y, i_t - 1));

if AHL(i_x, i_y, i_t) < 0
    AHL(i_x, i_y, i_t) = 0;
end

end
end

% ensure no flux out of plate
AHL(1, :, i_t) = AHL(2, :, i_t);
AHL(length(xgrid), :, i_t) = AHL(length(xgrid)-1, :, i_t);
AHL(:, 1, i_t) = AHL(:, 2, i_t);
AHL(:, length(ygrid), i_t) = AHL(:, length(xgrid)-1, i_t);

% update R, TX, and GFP

R1(:,:,i_t) = R1(:,:,i_t-1) + dt * (ror * highLux^2 .* AHL(:,:,i_t-1).^2 - DELr * R1(:,:,i_t-1));
R2(:,:,i_t) = R2(:,:,i_t-1) + dt * (ror * lowLux^2 .* AHL(:,:,i_t-1).^2 - DELr * R2(:,:,i_t-1));
R3(:,:,i_t) = R3(:,:,i_t-1) + dt * (ror * lowLux^2 .* AHL(:,:,i_t-1).^2 - DELr * R3(:,:,i_t-1));

TX1(:,:,i_t) = TX1(:,:,i_t-1) + dt * ((Atxgfp .* R1(:,:,i_t-1) / highKR)./(1+ R1(:,:,i_t-1) / highKR) - DELtxgfp .* TX1(:,:,i_t-1));
TX2(:,:,i_t) = TX2(:,:,i_t-1) + dt * ((Atxgfp .* R2(:,:,i_t-1) / highKR)./(1+ R2(:,:,i_t-1) / highKR) - DELtxgfp .* TX2(:,:,i_t-1));
TX3(:,:,i_t) = TX3(:,:,i_t-1) + dt * ((Atxgfp .* R3(:,:,i_t-1) / lowKR)./(1+ R3(:,:,i_t-1) / lowKR) - DELtxgfp .* TX3(:,:,i_t-1));

GFP1(:,:,i_t) = GFP1(:,:,i_t-1) + dt * (Agfp * TX1(:,:,i_t-1) - DELgfp * GFP1(:,:,i_t-1));
GFP2(:,:,i_t) = GFP2(:,:,i_t-1) + dt * (Agfp * TX2(:,:,i_t-1) - DELgfp * GFP2(:,:,i_t-1));
GFP3(:,:,i_t) = GFP3(:,:,i_t-1) + dt * (Agfp * TX3(:,:,i_t-1) - DELgfp * GFP3(:,:,i_t-1));

% get AHL and GFP edge

EDGEAHL1(i_t) = getSwitch(AHL(mid,1:mid,i_t), s1switch, xgrid);
EDGEAHL2(i_t) = getSwitch(AHL(mid,1:mid,i_t), s2switch, xgrid);
EDGEAHL3(i_t) = getSwitch(AHL(mid,1:mid,i_t), s3switch, xgrid);

EDGEGFP1(i_t) = getSwitch(GFP1(mid,1:mid,i_t), gfpswitch, xgrid);
EDGEGFP2(i_t) = getSwitch(GFP2(mid,1:mid,i_t), gfpswitch, xgrid);
EDGEGFP3(i_t) = getSwitch(GFP3(mid,1:mid,i_t), gfpswitch, xgrid);

end

%%

% finish up edge calculations

EDGEAHL1 = abs(EDGEAHL1);
EDGEAHL2 = abs(EDGEAHL2);
EDGEAHL3 = abs(EDGEAHL3);

EDGEGFP1 = abs(EDGEGFP1);
EDGEGFP2 = abs(EDGEGFP2);
EDGEGFP3 = abs(EDGEGFP3);

EDGEGFP1(EDGEGFP1 < 3) = 3;
EDGEGFP2(EDGEGFP2 < 3) = 3;
EDGEGFP3(EDGEGFP3 < 3) = 3;

% the GFP edge can't go backwards

[val, i] = max(EDGEGFP1);
EDGEGFP1(i:end) = val;

[val, i] = max(EDGEGFP2);
EDGEGFP2(i:end) = val;

[val, i] = max(EDGEGFP3);
EDGEGFP3(i:end) = val;

%% calculate NRMSE

data = readtable('gfp_data.xlsx');

% get model data from 0 - 24 hours

modS1 = EDGEGFP1(round((0:24) * 60 ./ dt)+1);
modS2 = EDGEGFP2(round((0:24) * 60 ./ dt)+1);
modS3 = EDGEGFP3(round((0:24) * 60 ./ dt)+1);

% get actual data and SEM

trueS1 = mean(table2array(data(:,2:4)), 2, 'omitnan');
trueS2 = mean(table2array(data(:,5:7)), 2, 'omitnan');
trueS3 = mean(table2array(data(:,8:10)), 2, 'omitnan');
SEMS1 = std(table2array(data(:,2:4)), 0,2, 'omitnan') / sqrt(3);
SEMS2 = std(table2array(data(:,5:7)), 0,2, 'omitnan') / sqrt(3);
SEMS3 = std(table2array(data(:,8:10)), 0,2, 'omitnan') / sqrt(3);

% calculate NRMSE

s1cost = NRMSE(modS1', trueS1);
s2cost = NRMSE(modS2', trueS2);
s3cost = NRMSE(modS3', trueS3);

% plot edge over time

figure(7)
tiledlayout(3,1)
nexttile
plot(0:24, trueS1, 'Color', 'b', 'LineWidth', 2)
hold on
plot(0:24, modS1, 'Color', 'r', 'LineWidth', 2)
plot(time./60, EDGEAHL1, 'LineWidth', 2)
errorbar(0:24, trueS1, SEMS1, 'Color', 'b', 'LineWidth', 1.5)
hold off
xlabel('Time (hours)','FontSize', 22)
ylabel('Distance (mm)','FontSize', 22)
xlim([0 24])
legend({'Data','Model', 'AHL'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 1 - NRMSE : %.3f', s1cost), 'FontSize', 22)

nexttile
plot(0:24, trueS2, 'Color', 'b', 'LineWidth', 2)
hold on
plot(0:24, modS2, 'Color', 'r', 'LineWidth', 2)
plot(time./60, EDGEAHL2, 'LineWidth', 2)
errorbar(0:24, trueS2, SEMS2, 'Color', 'b', 'LineWidth', 1.5)
hold off
xlabel('Time (hours)','FontSize', 22)
ylabel('Distance (mm)','FontSize', 22)
xlim([0 24])
legend({'Data','Model', 'AHL'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 2 - NRMSE : %.3f', s2cost), 'FontSize', 22)

nexttile
plot(0:24, trueS3, 'Color', 'b', 'LineWidth', 2)
hold on
plot(0:24, modS3, 'Color', 'r', 'LineWidth', 2)
plot(time./60, EDGEAHL3, 'LineWidth', 2)
errorbar(0:24, trueS3, SEMS3, 'Color', 'b', 'LineWidth', 1.5)
hold off
xlabel('Time (hours)','FontSize', 22)
ylabel('Distance (mm)','FontSize', 22)
xlim([0 24])
legend({'Data','Model', 'AHL'}, 'FontSize', 18, 'Location', 'northwest')
title(sprintf('Strain 3 - NRMSE : %.3f', s3cost), 'FontSize', 22)

%%

map = linspace(.1, 1, 100);
map = [zeros(100,1), map', zeros(100,1)];

% pull in images and find image center

img1 = imread('plate_images/strain1_rep1/1.0s_Exposure/2021_01_22_040011_1.0s_Exposure.png');
img2 = imread('plate_images/strain1_rep1/1.0s_Exposure/2021_01_22_120010_1.0s_Exposure.png');

[r1, c1, ppm1, ~] = find_center(rgb2gray(img1));
[r2, c2, ppm2, ~] = find_center(rgb2gray(img2));

% plot ahl concentration and gfp/image edges

figure(1)
subplot(2,3,1)
imagesc(xgrid, ygrid, AHL(:,:, round(8 * 60 / dt)))
ax = gca;
ylabel('T = 8 hours','FontSize', 24)
axis off
set(get(ax,'YLabel'),'Visible','on')
colorbar('FontSize', 14)
title('[AHL] uM','FontSize', 24)
ax1 = subplot(2,3,2);
imagesc(xgrid, ygrid, GFP1(:,:, round(8 * 60 / dt)))
axis off
colormap(ax1, map)
colorbar('FontSize', 14)
xline(EDGEGFP1(round(8 * 60 / dt)),'r', 'LineWidth', 2)
xline(-1*EDGEGFP1(round(8 * 60 / dt)),'r', 'LineWidth', 2)
title('[GFP] Model uM','FontSize', 24)
subplot(2,3,3)
imshow(img1((r1-500):(r1+500),(c1-500):(c1+500),:))
hold on
xline(500 + 256, 'r', 'LineWidth', 2)
xline(500 - 256, 'r', 'LineWidth', 2)
title('GFP Data','FontSize', 24)

subplot(2,3,4)
imagesc(xgrid, ygrid, AHL(:,:, round(16 * 60 / dt)))
ax = gca;
ylabel('T = 16 hours','FontSize', 24)
axis off
set(get(ax,'YLabel'),'Visible','on')
colorbar('FontSize', 14)
title('[AHL] uM','FontSize', 24)
ax2 = subplot(2,3,5);
imagesc(xgrid, ygrid, GFP1(:,:, round(16 * 60 / dt)))
axis off
xline(EDGEGFP1(round(16 * 60 / dt))-1,'r','LineWidth', 2)
xline(-1*EDGEGFP1(round(16 * 60 / dt))+1,'r', 'LineWidth', 2)
colorbar('FontSize', 14)
colormap(ax2, map)
title('[GFP] Model uM','FontSize', 24)
subplot(2,3,6)
imshow(img2((r2-300):(r2+700),(c2-290):(c2+690),:))
hold on
xline(500 + 339, 'r', 'LineWidth', 2)
xline(500 - 339, 'r', 'LineWidth', 2)
title('GFP Data','FontSize', 24)

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

function i_x = getSwitch(rows, switchConc, xgrid)

% find the first time that the concentration in rows is larger than the
% switching concentration

i = find(rows >= switchConc, 1, 'first');

if isempty(i) && (sum(rows) ~= 0)
    [~, i] = max(rows);
elseif isempty(i) && (sum(rows) == 0)
    i = length(rows);
end

i_x = xgrid(i);

end

function cost = NRMSE(Yexp, Y)

cost = sqrt(sum((Yexp - Y).^2)/length(Y)) / (max(Y) - min(Y));

end