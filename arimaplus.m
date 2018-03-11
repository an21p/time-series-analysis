clear all
close all

%% Input Data
% files = dir('data/output/*.csv')
% for f = files
%     fprintf("%s\n",f.name)
% end
cryptos_daily_close = getCrypo('data/output/cryptos_daily_close.csv', 'yyyy-MM-dd');
cryptos_daily_lr = getCrypo('data/output/cryptos_daily_lr.csv', 'yyyy-MM-dd');

cryptos_weekly_close = getCrypo('data/output/cryptos_weekly_close.csv', 'yyyy-MM-dd');
cryptos_weekly_lr = getCrypo('data/output/cryptos_weekly_lr.csv', 'yyyy-MM-dd');

cryptos_monthly_close = getCrypo('data/output/cryptos_monthly_close.csv', 'MMM yyyy');
cryptos_monthly_lr = getCrypo('data/output/cryptos_monthly_lr.csv', 'MMM yyyy');

fx_daily_close = getCrypo('data/output/fx_daily_close.csv', 'yyyy-MM-dd');
fx_daily_lr = getCrypo('data/output/fx_daily_lr.csv', 'yyyy-MM-dd');

fx_weekly_close = getCrypo('data/output/fx_weekly_close.csv', 'yyyy-MM-dd');
fx_weekly_lr = getCrypo('data/output/fx_weekly_lr.csv', 'yyyy-MM-dd');

fx_monthly_close = getCrypo('data/output/fx_monthly_close.csv', 'MMM yyyy');
fx_monthly_lr = getCrypo('data/output/fx_monthly_lr.csv', 'MMM yyyy');

%% Assets
% Cryptos
NEM_d = cryptos_daily_close.NEM;
NEM_dlc = log(NEM_d);
NEM_dlr = cryptos_daily_lr.NEM;

NEM_w = cryptos_weekly_close.NEM;
NEM_wlc = log(NEM_w);
NEM_wlr = cryptos_weekly_lr.NEM;

NEM_m = cryptos_monthly_close.NEM;
NEM_mlc = log(NEM_m);
NEM_mlr = cryptos_monthly_lr.NEM;


% FX time series
FX_d = fx_daily_close.Close;
FX_dlc = log(FX_d);
FX_dlr = fx_daily_lr.Log_Returns;

FX_w = fx_weekly_close.Close;
FX_wlc = log(FX_w);
FX_wlr = fx_weekly_lr.Log_Returns;

FX_m = fx_monthly_close.Close;
FX_mlc = log(FX_m);
FX_mlr = fx_monthly_lr.Log_Returns;


%% ARIMA Crypto

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,NEM_d);

% t-student Distribution
tdist = struct('Name','t','DoF',2.5); % DoF = 2.5 comes from part 1 of the assignment
Mdlt = arima(1,1,0);
Mdlt.Distribution = tdist;
EstMdlt = estimate(Mdlt,NEM_d);

% Forcasting
test = 250;
NEM_train = NEM_dlc(1:(end-test));
y = NEM_train;
T = length(y);
Mdlf = arima(1,1,0);
Mdlf.Distribution = tdist;
EstMdlf = estimate(Mdlf,y);

[yF,yMSE] = forecast(EstMdlf,test,'Y0',y);
upper = yF + 1.96*sqrt(yMSE);
lower = yF - 1.96*sqrt(yMSE);

figure
plot(NEM_dlc,'Color',[.75,.75,.75])
hold on
plot(NEM_train,'Color',[0,0,.75])
h1 = plot(T+1:T+test,yF,'r','LineWidth',2);
h2 = plot(T+1:T+test,upper,'k--','LineWidth',1.5);
plot(T+1:T+test,lower,'k--','LineWidth',1.5)
xlim([0,T+test])
title('Forecast and 95% Forecast Interval')
legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
hold off

%% CV
% 
% n = 10; % 200;
% pmseCR = zeros(n,1);
% for h = 1:n
%     fprintf("#")
%     y_train = NEM_dlc(1:(end-h));
%     y_test = NEM_dlc((end-h+1):end);
%     Mdlf = arima(1,1,0);
%     Mdlf.Distribution = tdist;
%     EstMdlf = estimate(Mdlf,y_train, 'Display', 'off');
%     [yF,~] = forecast(EstMdlf,h,'Y0',y_train);
%     pmseCR(h) = mean((y_test-yF).^2);
% end
% fprintf("\n")
% 
% %% CV Plot
% figure,
% plot(1:n, pmseCR, 'b--', 'LineWidth', 1.2)
% title('Prediction Error')


%% ARIMA FX

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,FX_d);

% t-student Distribution
tdist = struct('Name','t','DoF',5.5564); % DoF = 5.5564 comes from part 1 of the assignment
Mdlt = arima(1,1,0);
Mdlt.Distribution = tdist;
EstMdlt = estimate(Mdlt,FX_d);

% Forcasting
test = 500;
FX_train = FX_dlc(1:(end-test));
y = FX_train;
T = length(y);
Mdlf = arima(1,1,0);
Mdlf.Distribution = tdist;
EstMdlf = estimate(Mdlf,FX_train);

[yF,yMSE] = forecast(EstMdlf,test,'Y0',y);
upper = yF + 1.96*sqrt(yMSE);
lower = yF - 1.96*sqrt(yMSE);

figure
plot(FX_dlc,'Color',[.75,.75,.75])
hold on
plot(FX_train,'Color',[0,0,.75])
h1 = plot(T+1:T+test,yF,'r','LineWidth',2);
h2 = plot(T+1:T+test,upper,'k--','LineWidth',1.5);
plot(T+1:T+test,lower,'k--','LineWidth',1.5)
xlim([0,T+test])
title('Forecast and 95% Forecast Interval')
legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
hold off

%% CV
% 
% n = 10; % 200;
% pmseFX = zeros(n,1);
% for h = 1:n
%     fprintf("#")
%     y_train = FX_dlc(1:(end-h));
%     y_test = FX_dlc((end-h+1):end);
%     Mdlf = arima(1,1,0);
%     Mdlf.Distribution = tdist;
%     EstMdlf = estimate(Mdlf,y_train, 'Display', 'off');
%     [yF,~] = forecast(EstMdlf,h,'Y0',y_train);
%     pmseFX(h) = mean((y_test-yF).^2);
% end
% fprintf("\n")
% 
% %% CV Plot
% figure,
% plot(1:n, pmseFX, 'b--', 'LineWidth', 1.2)
% title('Prediction Error')

%% Stationarity Crypto

my_stationarity("NEM Daily", NEM_dlc, NEM_dlr); 
my_stationarity("NEM Weekly", NEM_wlc, NEM_wlr); 
my_stationarity("NEM Monthly", NEM_mlc, NEM_mlr); 

% ac_stationarity(NEM_dlc, NEM_wlc, NEM_mlc); 

%% Stationarity FX

my_stationarity("NEM Daily", FX_dlc, FX_dlr); 
my_stationarity("NEM Weekly", FX_wlc, FX_wlr); 
my_stationarity("NEM Monthly", FX_mlc, FX_mlr); 

% ac_stationarity(FX_dlc, FX_wlc, FX_mlc); 

%% Tail Crypto

neg_lr = sort(-NEM_dlr(NEM_dlr<0));
neg_rank = 1:length(neg_lr); % create rank of negative log returns
neg_ranked = 1 - (neg_rank/(length(neg_lr)+1));

tls = fitdist(NEM_dlr,'tlocationscale');
normdist = fitdist(NEM_dlr,'Normal');

m = normdist.mu;
s = normdist.sigma;

obs = 10;
len = size(neg_lr(:,1));
extreme = neg_lr((len-obs):len);
extreme_neg_rank = 1:length(extreme); % create rank of negative log returns
extreme_neg_ranked = 1 - (extreme_neg_rank/(length(extreme)+1));
ratio = extreme_neg_ranked/neg_ranked((len-obs):len);

f = fit(extreme,extreme_neg_ranked','b*x^(-alpha-1)', 'Start', [0 0])

xp = linspace(min(extreme),max(extreme)+0.2,100);
x = max(NEM_dlr)/1000:max(NEM_dlr)/1000:max(NEM_dlr);
alpha_NEM = f.alpha;


figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(extreme, neg_ranked((len-obs):len), 'ob')
loglog(xp,f(xp)/ratio,'-g', 'LineWidth', 2)
loglog(x,1-(cdf('tLocationScale',x,tls.mu,tls.sigma,tls.nu)-.5)*2,'-k','linewidth',2)
legend({'neg ret','extreme neg ret', 'power law', 'tLocationScale'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['power law fit for extreme negarive returns with \alpha = ' num2str(alpha_NEM)],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)

%% Tail FX

neg_lr = sort(-FX_dlr(FX_dlr<0));
neg_rank = 1:length(neg_lr); % create rank of negative log returns
neg_ranked = 1 - (neg_rank/(length(neg_lr)+1));

tls = fitdist(FX_dlr,'tlocationscale');
normdist = fitdist(FX_dlr,'Normal');

m = normdist.mu;
s = normdist.sigma;

obs = 35;
len = size(neg_lr(:,1));
extreme = neg_lr((len-obs):len);
extreme_neg_rank = 1:length(extreme); % create rank of negative log returns
extreme_neg_ranked = 1 - (extreme_neg_rank/(length(extreme)+1));
ratio = extreme_neg_ranked/neg_ranked((len-obs):len);

f = fit(extreme,extreme_neg_ranked','b*x^(-alpha-1)', 'Start', [0 0])

xp = linspace(min(extreme),max(extreme)+0.2,100);
x = max(FX_dlr)/1000:max(FX_dlr)/1000:max(FX_dlr);
alpha_FX = f.alpha;


figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(extreme, neg_ranked((len-obs):len), 'ob')
loglog(xp,f(xp)/ratio,'-g', 'LineWidth', 2)
loglog(x,1-(cdf('tLocationScale',x,tls.mu,tls.sigma,tls.nu)-.5)*2,'-k','linewidth',2)
legend({'neg ret','extreme neg ret', 'power law', 'tLocationScale'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['power law fit for extreme negarive returns with \alpha = ' num2str(alpha_FX)],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
%% Hurst Exponent Crypto

% figure
% genhurst2(NEM_d,1:.2:3);
% 
% figure
% genhurst2(NEM_w,1:.2:3);

fprintf("\n\n");
fprintf("Hurst Exponend for %s: %f\n", "NEM Daily", genhurst(NEM_d,2));
fprintf("Hurst Exponend for %s: %f\n", "NEM Weekly", genhurst(NEM_w,2));

k_tau_plot("NEM", NEM_d, 10);

%% Hurst Exponent FX

% figure
% genhurst2(FX_d,1:.2:3);
% 
% figure
% genhurst2(FX_w,1:.2:3);

fprintf("\n\n");
fprintf("Hurst Exponend for %s: %f\n", "NEM Daily", genhurst(FX_d,2));
fprintf("Hurst Exponend for %s: %f\n", "NEM Weekly", genhurst(FX_w,2));

k_tau_plot("FX", FX_d, 10);

%% SCALING FROM ALPHA (DISTRIBUTION) Crypto % From Antonio
% from daily to weekly 
freq = 7;
hurst2 = genhurst(NEM_d,2)
hurst2 = 1/alpha_NEM
new_lret = NEM_dlr * freq^hurst2;


rescaled_pd = fitdist(new_lret,'Normal');
x_values = (-10:1:10);
y = pdf(rescaled_pd,x_values);
y = y*freq^hurst2;

pd = fitdist(NEM_wlr,'Normal');
theoretical = pdf(pd,x_values);

figure
hold on
plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
plot(x_values,log(theoretical/max(theoretical)), 'x--', 'LineWidth',2)
legend('rescaled from daily','weekly theoretical', 'location', 'best')
hold off
KLDiv(theoretical,y)

% from daily to monthly 
freq = 30;
new_lret = NEM_dlr * freq^hurst2;


rescaled_pd = fitdist(new_lret,'Normal');
x_values = (-10:1:10);
y = pdf(rescaled_pd,x_values);
y = y*freq^hurst2;

pd = fitdist(NEM_mlr,'Normal');
theoretical = pdf(pd,x_values);

figure
hold on
plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
plot(x_values,log(theoretical/max(theoretical)), 'x--', 'LineWidth',2)
legend('rescaled from daily','monthly theoretical', 'location', 'best')
hold off
KLDiv(theoretical,y)

%% SCALING FROM ALPHA (DISTRIBUTION) FX % From Antonio
% from daily to weekly 
freq = 7;
hurst2 = genhurst(FX_d,2)
hurst2 = 1/alpha_FX
new_lret = FX_dlr * freq^hurst2;


rescaled_pd = fitdist(new_lret,'Normal');
x_values = (-10:.1:10);
y = pdf(rescaled_pd,x_values);
y = y*freq^hurst2;

pd = fitdist(FX_wlr,'Normal');
theoretical = pdf(pd,x_values);

figure
hold on
plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
plot(x_values,log(theoretical/max(theoretical)), 'x--', 'LineWidth',2)
legend('rescaled from daily','weekly theoretical', 'location', 'best')
hold off
KLDiv(theoretical,y)

% from daily to monthly 
freq = 30;
new_lret = FX_dlr * freq^hurst2;

rescaled_pd = fitdist(new_lret,'Normal');
x_values = (-10:.1:10);
y = pdf(rescaled_pd,x_values);
y = y*freq^hurst2;

pd = fitdist(FX_mlr,'Normal');
theoretical = pdf(pd,x_values);

figure
hold on
plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
plot(x_values,log(theoretical/max(theoretical)), 'o--', 'LineWidth',2)
legend('rescaled from daily','monthly theoretical', 'location', 'best')
hold off
KLDiv(theoretical,y)
