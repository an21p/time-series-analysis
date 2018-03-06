clear all
close all

%% Input Data
% files = dir('data/output/*.csv')
% for f = files
%     fprintf("%s\n",f.name)
% end
cryptos_daily_log_close = getCrypo('data/output/cryptos_daily_log_close.csv', 'yyyy-MM-dd');
cryptos_daily_lr = getCrypo('data/output/cryptos_daily_lr.csv', 'yyyy-MM-dd');

cryptos_weekly_log_close = getCrypo('data/output/cryptos_weekly_log_close.csv', 'yyyy-MM-dd');
cryptos_weekly_lr = getCrypo('data/output/cryptos_weekly_lr.csv', 'yyyy-MM-dd');

cryptos_monthly_log_close = getCrypo('data/output/cryptos_monthly_log_close.csv', 'MMM yyyy');
cryptos_monthly_lr = getCrypo('data/output/cryptos_monthly_lr.csv', 'MMM yyyy');


%% ARIMA 

NEM = cryptos_daily_log_close.NEM;

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,NEM);

% t-student Distribution
tdist = struct('Name','t','DoF',2.5); % DoF = 2.5 comes from part 1 of the assignment
Mdlt = arima(1,1,0);
Mdlt.Distribution = tdist;
EstMdlt = estimate(Mdlt,NEM);


% Forcasting
test = 200;
NEM_train = NEM(1:(end-test));
y = NEM_train;
T = length(y);
Mdlf = arima(1,1,0);
Mdlf.Distribution = tdist;
EstMdlf = estimate(Mdlf,NEM_train);

[yF,yMSE] = forecast(EstMdlf,test,'Y0',y);
upper = yF + 1.96*sqrt(yMSE);
lower = yF - 1.96*sqrt(yMSE);

figure
hold on
plot(NEM,'Color',[.75,.75,.75])
h1 = plot(T+1:T+test,yF,'r','LineWidth',2);
h2 = plot(T+1:T+test,upper,'k--','LineWidth',1.5);
plot(T+1:T+test,lower,'k--','LineWidth',1.5)
xlim([0,T+test])
title('Forecast and 95% Forecast Interval')
legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
hold off


%% Functions

function t = getCrypo(filename, dateformat)
    table = readtable(filename);
    table.Var1 = datetime(table.Var1,'InputFormat',dateformat);
    table.Properties.VariableNames(1) = {'Date'};
    t = table;
end
