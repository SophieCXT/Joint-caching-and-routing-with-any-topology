% Online setting tests

%Import data
t = csvread('data/weeks_views.csv',1,0);% = readtable('data/weeks_views.csv');
%weeks_views = table2dataset(t);
window_size = 2;
for video = 1:50%length(weeks_views), top 50 (C)
    views_series = t(video,:).';
    exp_ts(video,:) = tsmovavg(views_series,'e',window_size,1);
%      triangular moving average moving average.
end


cmap = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560;...
    0.4660, 0.6740, 0.1880;...
    0.3010, 0.7450, 0.9330;...
    0.6350, 0.0780, 0.1840];
markers = {'.','.','.','.','.','.','.','.' };%{'p','o','s','d','x','^','+'};
size1=16; % label font size
size2=14; % legend font size
step = 3; % step size in markers
markersize = 11;

figure;
hold on;
for i=1:7
    plot(t(i,window_size:20),'Color', cmap(i,:),'LineWidth', 1.5);
end
for i=1:7
    plot(exp_ts(i,window_size:20),'Color', cmap(i,:),'LineWidth', 1,'LineStyle','--');
end
xlabel('weeks', 'FontSize', size1)
ylabel('#views', 'FontSize', size1)
title(['EWMA Prediction VS Original Time Series, where n=' num2str(window_size)])

%% old model
% % Let t be a 1-D array of data points
% lambda = 0.2; %forgetting factor, n=9,7,5,3
% prevWeightFactor = 0; %initialise weight factor and moving average
% prevAvg = 0;
% 
% for video = 1:12130%length(t)
% for i = 1:12 % test set 60% 
%     presentWeightFactor = lambda * prevWeightFactor + 1;
%     presentAvg = (1 - (1/presentWeightFactor)) * prevAvg + (1/presentWeightFactor) * t(i);
%     if (presentAvg < 0.5 * prevAvg) || (presentAvg > 1.5 * prevAvg)
%         presentAvg = prevAvg;  %ignore this input, you might want to skip this step for the first sample
%     else   %accept this input in the moving average
%         prevWeightFactor = presentWeightFactor;
%         prevAvg = presentAvg;
%     end
% end
% end
