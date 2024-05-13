close all;
clear;
h = [1.2 1.5 2.0 5.5];
sigma = [0.05 0.1 0.25 0.1]*2;

longitude = [45.50 45.25 45.02 44.90];

% Specify the file name and path
% filename = "mcmc_Kane_mantle";
filename = "mcmc_Kane_mantle_tg_hmin_30.000000_new";

h_min = 30;
H = 60;

% Open the file for reading
fileID = fopen(filename, 'r');

% Read the first four rows
header = textscan(fileID, '%f', 4);
values = header{1};
numRows = values(1);
numColsA = values(2);
numColsB = values(3);
numColsC = values(4);

% Read the remaining data
data = textscan(fileID, repmat('%f', 1, numColsA+numColsB+numColsC), 'Delimiter', ' ');

% Close the file
fclose(fileID);

% Extract the A, B, and C matrices from the data
var_list = cell2mat(data(:, 1:numColsA));
y_list = cell2mat(data(:, numColsA+1:numColsA+numColsB));
cost_list = cell2mat(data(:, numColsA+numColsB+1:end));

% index_select = cost_list<0.5;

FAD = var_list(:, 1);
XRP_max = var_list(:, 2);
xrp_drops = var_list(:, 3:5);
XRP_segments = var_list(:, 6:9);
F_max = var_list(:, 10);

segments_length_total = 60;
total_segments = sum(XRP_segments, 2);
XRP_segments = XRP_segments ./ total_segments * segments_length_total;
segments = [0*XRP_segments(:, 1), cumsum(XRP_segments(:, 1:end-1), 2)];

total_xrp_drops = sum(xrp_drops, 2);
xrp_drops = xrp_drops ./ total_xrp_drops .* XRP_max;
xrps = zeros(numRows, 4);
xrps(:, 1) = XRP_max;
for ii = 2:4
    xrps(:, ii) = xrps(:, ii-1) - xrp_drops(:, ii-1);
end

segments = [segments(:, 1)*0, segments+60, segments(:, 1)*0+140];
xrps(xrps<0) = 0;
xrps = [xrps(:,1), xrps, xrps(:, 1)*0];
%%
hplot = figure(1);
set(hplot,'Position',[20,200,300,300]);

plot(-longitude, y_list, 'r:'); hold on;
errorbar(-longitude, h, sigma, 'kd');
title('Crustal thickness (km)');

hplot2 = figure(2);
set(hplot2,'Position',[300,200,900,600]);
subplot(2,3,[1, 2, 3]);
for ii=1:numRows
    plot(segments(ii, :), xrps(ii, :), 'r:'); hold on;
end

xlabel('Distance away from 45.5W (km)');
title('Refractory peridotite vol.% from 45.5W to axis');

subplot(2,3,4);
data_hist = FAD;
mean(data_hist)
std(data_hist)
histogram(data_hist); hold on;
title('Degree of ancient depletion of RP');


subplot(2,3,5);
data_hist = F_max*(H-h_min)/H - FAD;
mean(data_hist)
std(data_hist)
histogram(data_hist); hold on;
title('Degree of melting of RP');

subplot(2,3,6);
data_hist = F_max*(H-h_min)/H;
mean(data_hist)
std(data_hist)
histogram(data_hist);
title('Degree of melting of FP');

hplot3 = figure(3);
set(hplot3,'Position',[200,400,300,300]);
plot(cost_list, 'k-'); hold on;
title('Misfit convergence');
%%
hplot4 = figure(4);
set(hplot4,'Position',[800,200,800,800]);
subplot(1,3,1);
for ii=1:numRows
    plot(xrps(ii, :), segments(ii, :), 'r:'); hold on;
end
set(gca, 'ydir','reverse');
xlabel('Refractory peridotite vol.%');

subplot(1,3,[2, 3]);
radius_earth = 6400;
longitude = [45.50 45.25 45.02 44.90];
ridge_azimuthal_angle = 5/180*pi; % NNE 5 degree

for ii=1:length(longitude)
    long = longitude(ii);
    z = (longitude(ii) - longitude(1))*radius_earth*pi/180*cos(ridge_azimuthal_angle);
    plot([0, 60, -60, 0], [0, 60, 60, 0]+z, 'b-'); hold on;
    plot([30, -30], [30, 30]+z, 'b--'); hold on;
end
set(gca, 'ydir','reverse');
ylim([0, 140]);

