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
Mdln = arima(2,1,0);
EstMdln = estimate(Mdln,NEM_dlc);

% t-student Distribution
tls = fitdist(NEM_dlr,'tlocationscale');
tdist = struct('Name','t','DoF',tls.nu);

% Forcasting
test = 259;
NEM_train = NEM_dlc(1:(end-test));
y = NEM_train;
T = length(y);
Mdlf = arima(2,1,0);
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


%% ARIMA FX daily

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,FX_dlc);

% t-student Distribution
tls = fitdist(FX_dlr,'tlocationscale');
tdist = struct('Name','t','DoF',tls.nu);

% Forcasting
test = 839;
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

%% ARIMA FX weekly

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,FX_wlc);

% t-student Distribution
tls = fitdist(FX_wlr,'tlocationscale');
tdist = struct('Name','t','DoF',tls.nu);


% Forcasting
test = 100;
FX_train = FX_wlc(1:(end-test));
y = FX_train;
T = length(y);
Mdlf = arima(1,1,0);
Mdlf.Distribution = tdist;
EstMdlf = estimate(Mdlf,FX_train);

[yF,yMSE] = forecast(EstMdlf,test,'Y0',y);
upper = yF + 1.96*sqrt(yMSE);
lower = yF - 1.96*sqrt(yMSE);

figure
plot(FX_wlc,'Color',[.75,.75,.75])
hold on
plot(FX_train,'Color',[0,0,.75])
h1 = plot(T+1:T+test,yF,'r','LineWidth',2);
h2 = plot(T+1:T+test,upper,'k--','LineWidth',1.5);
plot(T+1:T+test,lower,'k--','LineWidth',1.5)
xlim([0,T+test])
title('Forecast and 95% Forecast Interval')
legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
hold off


%% ARIMA FX monthly

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,FX_mlc);

% t-student Distribution
tls = fitdist(FX_mlr,'tlocationscale');
tdist = struct('Name','t','DoF',tls.nu);


% Forcasting
test = 20;
FX_train = FX_mlc(1:(end-test));
y = FX_train;
T = length(y);
Mdlf = arima(1,1,0);
Mdlf.Distribution = tdist;
EstMdlf = estimate(Mdlf,FX_train);

[yF,yMSE] = forecast(EstMdlf,test,'Y0',y);
upper = yF + 1.96*sqrt(yMSE);
lower = yF - 1.96*sqrt(yMSE);

figure
plot(FX_mlc,'Color',[.75,.75,.75])
hold on
plot(FX_train,'Color',[0,0,.75])
h1 = plot(T+1:T+test,yF,'r','LineWidth',2);
h2 = plot(T+1:T+test,upper,'k--','LineWidth',1.5);
plot(T+1:T+test,lower,'k--','LineWidth',1.5)
xlim([0,T+test])
title('Forecast and 95% Forecast Interval')
legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
hold off


%% ARIMA FX daily ?exp()?

% Normal Distribution
Mdln = arima(1,1,0);
EstMdln = estimate(Mdln,FX_dlc);

% t-student Distribution
tls = fitdist(FX_dlr,'tlocationscale');
tdist = struct('Name','t','DoF',tls.nu);

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
yF = exp(yF);
upper = exp(upper);
lower = exp(lower);

figure
plot(FX_d,'Color',[.75,.75,.75])
hold on
plot(exp(FX_train),'Color',[0,0,.75])
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


neg_lr_w = sort(-NEM_wlr(NEM_wlr<0)); %/7^genhurst(NEM_d,2);
neg_rank_w = 1:length(neg_lr_w); % create rank of negative log returns
neg_ranked_w = 1 - (neg_rank_w/(length(neg_lr_w)+1));

neg_lr_m = sort(-NEM_mlr(NEM_mlr<0)); %/30^genhurst(NEM_d,2);
neg_rank_m = 1:length(neg_lr_m); % create rank of negative log returns
neg_ranked_m = 1 - (neg_rank_m/(length(neg_lr_m)+1));

tls = fitdist(NEM_dlr,'tlocationscale');
normdist = fitdist(NEM_dlr,'Normal');

m = normdist.mu;
s = normdist.sigma;

obs = 11;
len = size(neg_lr(:,1));
extreme = neg_lr((len-obs):len);
extreme_neg_rank = 1:length(extreme); % create rank of negative log returns
extreme_neg_ranked = 1 - (extreme_neg_rank/(length(extreme)+1));
ratio = extreme_neg_ranked/neg_ranked((len-obs):len);

f = fit(extreme,extreme_neg_ranked','b*x^(-alpha-1)', 'Start', [0 0])
alpha_NEM = f.alpha;

xp = linspace(min(extreme),max(extreme)+0.2,100);
x = max(NEM_dlr)/1000:max(NEM_dlr)/1000:max(NEM_dlr);

figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(extreme, neg_ranked((len-obs):len), 'ob')
loglog(xp,f(xp)/ratio,'-g', 'LineWidth', 2)
loglog(x,1-(cdf('tLocationScale',x,tls.mu,tls.sigma,tls.nu)-.5)*2,'-k','linewidth',2)
legend({'daily neg ret','extreme neg ret', 'power law', 'tLocationScale'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['power law fit for extreme negarive returns with \alpha = ' num2str(alpha_NEM)],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)

figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(neg_lr_w, neg_ranked_w, 'xb')
loglog(neg_lr_m, neg_ranked_m, 'xg')
legend({'daily neg ret','weekly neg ret','monthly neg ret'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['NEM Log-Returns at different time-horizons'],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)

%% Tail FX

neg_lr = sort(-FX_dlr(FX_dlr<0));
neg_rank = 1:length(neg_lr); % create rank of negative log returns
neg_ranked = 1 - (neg_rank/(length(neg_lr)+1));

neg_lr_w = sort(-FX_wlr(FX_wlr<0)); %/7^genhurst(FX_d,2);
neg_rank_w = 1:length(neg_lr_w); % create rank of negative log returns
neg_ranked_w = 1 - (neg_rank_w/(length(neg_lr_w)+1));

neg_lr_m = sort(-FX_mlr(FX_mlr<0)); %/30^genhurst(FX_d,2);;
neg_rank_m = 1:length(neg_lr_m); % create rank of negative log returns
neg_ranked_m = 1 - (neg_rank_m/(length(neg_lr_m)+1));

tls = fitdist(FX_dlr,'tlocationscale');
normdist = fitdist(FX_dlr,'Normal');

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
alpha_FX = f.alpha;

xp = linspace(min(extreme),max(extreme)+0.2,100);
x = max(FX_dlr)/1000:max(FX_dlr)/1000:max(FX_dlr);

figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(extreme, neg_ranked((len-obs):len), 'ob')
loglog(xp,f(xp)/ratio,'-g', 'LineWidth', 2)
loglog(x,1-(cdf('tLocationScale',x,tls.mu,tls.sigma,tls.nu)-.5)*2,'-k','linewidth',2)
legend({'daily neg ret','extreme neg ret', 'power law', 'tLocationScale'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['power law fit for extreme negarive returns with \alpha = ' num2str(alpha_FX)],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)

figure,
loglog(neg_lr, neg_ranked, 'xr')
hold on
loglog(neg_lr_w, neg_ranked_w, 'xb')
loglog(neg_lr_m, neg_ranked_m, 'xg')
legend({'daily neg ret','weekly neg ret','monthly neg ret'},'Location','best')
axis([1e-5 1 1e-4 1])
title(['FX Log-Returns at different time-horizons'],'fontsize',14)
xlabel('log-return','fontsize',14)
ylabel('complemetary cumulative distribution','fontsize',14)
set(gca,'fontsize',14)
%% Hurst Exponent Crypto

figure
genhurst2(NEM_d,1:.1:3);
xlim([1 20])

fprintf("\n\n");
fprintf("Hurst Exponend for %s: %f\n", "NEM Daily", genhurst(NEM_d,2));

k_tau_plot("NEM", NEM_d, 10);

%% Hurst Plot on P.73
S=NEM_d;
q = [1:0.01:4];
maxT=19;
h_q = genhurst(S,q,maxT);
qhq = h_q.*q';
figure,
plot(q,qhq)
hold on
title('Hurst Exponents with Moments'); 
xlabel('q');
ylabel('qH(q)');

%% Stylized Fact 6 
pr=NEM_d;
syb = 'v<>^';
col = 'rbmg';
k=0;
mm=[1 7 30];
figure
for m=mm
     k=k+1;
     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
     [f,b]=hist(lr,20);
     %subplot(1,2,1)
     semilogy(b,f/max(f),['-',syb(k),col(k)])
     hold on
end
 axis([-1.5 +1.5 1e-4 1.1])
 legend(num2str(mm'),'Location','best')
 set(gca,'fontsize',12)
 xlabel('log-return','fontsize',12)
 ylabel('relative freq','fontsize',12)
 title('Unscaled of NEM Log Return Distribution')
 print('1secRetDistr.eps','-depsc')

%% 
 hh=genhurst(pr,2);
 k=0;
 figure
 for m=mm
     k=k+1;
     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
     %subplot(1,2,2)
     [f,b]=hist(lr,20);
     semilogy(b/m.^hh,f/max(f),['-',syb(k),col(k)])
     hold on
 end
 text(-1.05,0.05,['H(2) = ',num2str(hh)],'fontsize',12)
 axis([-1.5 +1.5 1e-4 1.1])
 legend(num2str(mm'),'Location','best')
 set(gca,'fontsize',12)
 xlabel('log-return','fontsize',12)
 ylabel('relative freq','fontsize',12)
 title('Rescaled NEM Log Return Distribution')
 print('1secRetDistrScaled.eps','-depsc')

%% Hurst Exponent FX

figure
genhurst2(FX_d,1:.1:3); 
xlim([1 20])

fprintf("\n\n");
fprintf("Hurst Exponend for %s: %f\n", "NEM Daily", genhurst(FX_d,2));

k_tau_plot("FX", FX_d, 10);


%% Hurst Plot on P.73
S=FX_d;
q = [1:0.01:4];
maxT=19;
h_q = genhurst(S,q,maxT);
qhq = h_q.*q';
figure,
plot(q,qhq)
title('Hurst Exponents with Moments'); 
xlabel('q');
ylabel('qH(q)');

%% Stylized Fact 6 
t1=find(diff(log(FX_d))~=0,1,'first')+1;
t2=find(diff(log(FX_d))~=0,1,'last');
pr=FX_d(t1:t2);
syb = 'v<>^';
col = 'rbmg';
k=0;
mm=[1 7 30];
figure
for m=mm
     k=k+1;
     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
     [f,b]=hist(lr,20);
     %subplot(1,2,1)
     semilogy(b,f/max(f),['-',syb(k),col(k)])
     hold on
 end
 axis([-0.2 +0.2 1e-4 1.1])
 legend(num2str(mm'),'Location','best')
 set(gca,'fontsize',12)
 xlabel('log-return','fontsize',12)
 ylabel('relative freq','fontsize',12)
 title('Unscaled of FX Log Return Distribution')
 print('1secRetDistr.eps','-depsc')

%% 
 hh=genhurst(pr,2);
 k=0;
 figure
 for m=mm
     k=k+1;
     lr = log(pr((m+1):m:end))-log(pr(1:m:(end-m)));
     %subplot(1,2,2)
     [f,b]=hist(lr,20);
     semilogy(b/m.^hh,f/max(f),['-',syb(k),col(k)])
     hold on
 end
 text(0.1,0.05,['H(2) = ',num2str(hh)],'fontsize',12)
 axis([-0.2 +0.2 1e-4 1.1])
 legend(num2str(mm'),'Location','best')
 set(gca,'fontsize',12)
 xlabel('log-return','fontsize',12)
 ylabel('relative freq','fontsize',12)
 title('Rescaled FX Log Return Distribution')
 print('1secRetDistrScaled.eps','-depsc')


%% SCALING FROM ALPHA (DISTRIBUTION) Crypto % From Antonio
% % % multiscaling 
% % % from weekly to daily 
% % freq = 7;
% % ghurst2 = genhurst(NEM_d,2)
% % alpha_NEM
% % hurst2 = 1/alpha_NEM
% % new_lret = NEM_wlr * (1/freq)^hurst2;
% % 
% % 
% % rescaled_pd = fitdist(new_lret,'Normal');
% % x_values = (-10:1:10);
% % y = pdf(rescaled_pd,x_values);
% % y = y*(1/freq)^hurst2;
% % 
% % pd = fitdist(NEM_dlr,'Normal');
% % theoretical = pdf(pd,x_values); 
% % 
% % figure
% % hold on
% % plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
% % plot(x_values,log(theoretical/max(theoretical)), 'x--', 'LineWidth',2)
% % legend('rescaled from weekly','daily theoretical', 'location', 'best')
% % xlim([-15 15])
% % hold off
% % KLDiv(theoretical,y)
% % 
% % hurst2 = 1/1.75;
% % ghurst1 = genhurst(NEM_d,1)
% % 
% % % from monthly to daily
% % freq = 30;
% % new_lret = NEM_mlr * (1/freq)^ghurst1;
% % 
% % rescaled_pd = fitdist(new_lret,'Normal');
% % x_values = (-10:1:10);
% % y = pdf(rescaled_pd,x_values);
% % y = y*(1/freq)^ghurst1;
% % 
% % pd = fitdist(NEM_dlr,'Normal');
% % theoretical = pdf(pd,x_values);
% % 
% % figure
% % hold on
% % plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
% % plot(x_values,log(theoretical/max(theoretical)), 'o--', 'LineWidth',2)
% % legend('rescaled from monthly','daily theoretical', 'location', 'best')
% % xlim([-15 15])
% % hold off
% % KLDiv(theoretical,y)

%% Scaling 

NEM_daily = fitdist(NEM_dlr,'Kernel');
NEM_weekly = fitdist(NEM_wlr,'Kernel');
NEM_monthly = fitdist(NEM_mlr,'Kernel');

x_values = linspace(-2,2,1000);
figure
h1 = plot(x_values, pdf(NEM_daily,x_values),'LineWidth',1.5);
hold on
h2 = plot(x_values, pdf(NEM_weekly,x_values),'LineWidth',1.5);
h3 = plot(x_values, pdf(NEM_monthly,x_values),'LineWidth',1.5);

ylabel('Probability Density')
xlabel('Log-Returns')
title(['NEM',' ','Log-Return Distribution'])
legend([h1,h2,h3],'Daily','Weekly','Monthly','Location','best')

%% SCALING FROM ALPHA (DISTRIBUTION) FX % From Antonio
% % % uniscaling 
% % % from weekly to daily 
% % freq = 7;
% % ghurst2 = genhurst(FX_d,2)
% % alpha_FX
% % hurst2 = 1/alpha_FX
% % new_lret = FX_wlr * (1/freq)^hurst2;
% % 
% % 
% % rescaled_pd = fitdist(new_lret,'Normal');
% % x_values = (-10:.1:10);
% % y = pdf(rescaled_pd,x_values);
% % y = y*(1/freq)^hurst2;
% % 
% % pd = fitdist(FX_dlr,'Normal');
% % theoretical = pdf(pd,x_values);
% % 
% % figure
% % hold on
% % plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
% % plot(x_values,log(theoretical/max(theoretical)), 'x--', 'LineWidth',2)
% % legend('rescaled from weekly','daily theoretical', 'location', 'best')
% % xlim([-1 1])
% % hold off
% % KLDiv(theoretical,y)
% % 
% % % from monthly to daily
% % freq = 30;
% % new_lret = FX_mlr * (1/freq)^hurst2;
% % 
% % rescaled_pd = fitdist(new_lret,'Normal');
% % x_values = (-10:.1:10);
% % y = pdf(rescaled_pd,x_values);
% % y = y*(1/freq)^hurst2;
% % 
% % pd = fitdist(FX_dlr,'Normal');
% % theoretical = pdf(pd,x_values);
% % 
% % figure
% % hold on
% % plot(x_values,log(y/max(y)), 'o--', 'LineWidth',2);
% % plot(x_values,log(theoretical/max(theoretical)), 'o--', 'LineWidth',2)
% % legend('rescaled from monthly','daily theoretical', 'location', 'best')
% % xlim([-1 1])
% % hold off
% % KLDiv(theoretical,y)

%% Scaling 

FX_daily = fitdist(FX_dlr,'Kernel');
FX_weekly = fitdist(FX_wlr,'Kernel');
FX_monthly = fitdist(FX_mlr,'Kernel');

x_values = linspace(-.5,.5,1000);
figure
h1 = plot(x_values, pdf(FX_daily,x_values),'LineWidth',1.5);
hold on
h2 = plot(x_values, pdf(FX_weekly,x_values),'LineWidth',1.5);
h3 = plot(x_values, pdf(FX_monthly,x_values),'LineWidth',1.5);
xlim([-.2 .2])
ylabel('Probability Density')
xlabel('Log-Returns')
title(['FX',' ','Log-Return Distribution'])
legend([h1,h2,h3],'Daily','Weekly','Monthly','Location','best')


%% Correlations
clear all
close all

cryptos = {'Bitcoin', 'Dash', 'Ethereum', 'Litecoin', 'Monero', 'NEM', 'Ripple', 'Siacoin', 'Stellar', 'Verge'};

all_crypto_table = readtable('data/output/all_crypto_log_returns.csv');
all_crypto = all_crypto_table{:,2:end};

PR = false;
KR = false;
SR = false;

% Pearson
title = "Pearson's Correlation Coefficients for 10 Cryptocurrencies (Daily Returns)";
if PR == false, [PR,~] = corr(all_crypto); end
colorcorr(title,PR,cryptos)

% Kendall
title = "Kendall's Rank Correlation Coefficients (\tau) for 10 Cryptocurrencies (Daily Returns)";
if KR == false, [KR,~] = corr(all_crypto,'type','Kendall'); end
colorcorr(title,KR,cryptos)

% Spearman
title = "Spearman?s Rank Correlation Coefficients (\rho) for 10 Cryptocurrencies (Daily Returns)";
if SR == false, [SR,~] = corr(all_crypto,'type','Spearman'); end
colorcorr(title,SR,cryptos)


%% PCA

[coeff,score,latent,tsquared,explained] = pca(zscore(all_crypto));
explained
cumsum(explained)

% figure
% scatter3(score(:,1),score(:,2),score(:,3))
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')

figure
subplot(1,2,1)
plot(explained, 'o--')
xlabel("Principal Component")
ylabel("Proportion of Variance Explained")

subplot(1,2,2)
plot(cumsum(explained), 'o--')
xlabel("Principal Component")
ylabel("Cumulative Proportion of Variance Explained")
ylim([0 100])

figure
biplot(coeff(:,[2 3]),'scores',score(:,[2 3]),'varlabels',cryptos);

% explained =
%    30.4380
%    12.1630
%     9.6110
%     8.8008
%     8.5273
%     7.5565
%     6.9467
%     6.4336
%     4.8366
%     4.6865

%% Hierarchical Clustering

[PR,~] = corr(zscore(all_crypto));
RPR = 1-PR;
E = pdist(RPR,'euclidean');
Y = squareform(E);
%Y = sqrt(2*(1-PR));

Z = linkage(Y,'complete','euclidean');

figure
dendrogram(Z,'ColorThreshold','default', 'Orientation', 'left' )
ind = str2num(get(gca,'YTickLabel'));
set(gca, 'YTickLabel',cryptos(ind))
ax = gca;
ax.Title.String = 'Hierarchical Clustering of Cryptocurrencies';
xtickangle(gca, 90)

c = cluster(Z,'maxclust',5);
figure
scatter3(Y(:,1),Y(:,2),Y(:,3),10,c, 'LineWidth', 5)
figure
scatter(Y(:,1),Y(:,2),10,c, 'LineWidth', 5)



%% MI

% Bitcoin = all_crypto_table.Bitcoin;
% Stellar = all_crypto_table.Stellar;
% Litecoin = all_crypto_table.Litecoin;
% Siacon = all_crypto_table.Siacon;
% NEM = all_crypto_table.NEM;
% 
% Ripple = all_crypto_table.Ripple;
% 
% 
% mi(Bitcoin,NEM)  % 0.0225
% mi(Stellar,NEM)  % 0.0339
% mi(Litecoin,NEM) % 0.0122
% mi(Siacon,NEM)   % 0.0392
% 
% mi(Litecoin, Bitcoin) % 0.1904
% mi(Ripple, Stellar) % 0.0726

%[F, cv] = granger_cause(Bitcoin, NEM, 0.05, 50)

%% Permutation Test

pvals = zeros(length(all_crypto(1,:)));
perms = 5000;

for i = 1:length(all_crypto(1,:))
    for j = 1:length(all_crypto(1,:))
        display(sprintf('i = %d, j = %d', i,j))
        pvals(i,j)=permtest(all_crypto(:,i), all_crypto(:,j), perms);
    end
end

%%

colorcorr("P values for Spearman Correlations from Permutation Test", pvals, cryptos, flipud(summer), 0,1)

%% Eigenvectors
