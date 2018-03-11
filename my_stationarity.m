function [] = my_stationarity(name, logPrices, logReturns)

fprintf("\n\nStationarity\n");
fprintf("--------------------------------\n");
% adftest    0 = unit root process | 1 = autoregressive
fprintf("adftest    0 = unit root process | 1 = autoregressive\n")
% kpsstest   0 = trend stationary | 1 = nonstationary unit root process
fprintf("kpsstest   0 = trend stationary  | 1 = nonstationary unit root process\n")
% vratiotest 0 = random walk | 1 = not random walk
fprintf("vratiotest 0 = random walk       | 1 = not random walk\n")
fprintf("--------------------------------\n");
fprintf("\n%s\n",name);
fprintf("----------------\n");
fprintf("\nLog Price\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(logPrices));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(logPrices));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(logPrices));
fprintf("\nLog Returns\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(logReturns));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(logReturns));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(logReturns));
end

