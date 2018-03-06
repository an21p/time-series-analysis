function [] = stationarity(name, d , w, m , dlr, wlr, mlr)

fprintf("\n\nStationarity\n");
fprintf("--------------------------------\n");
% adftest    0 = unit root process | 1 = autoregressive
fprintf("adftest    0 = unit root process | 1 = autoregressive\n")
% kpsstest   0 = trend stationary | 1 = nonstationary unit root process
fprintf("kpsstest   0 = trend stationary  | 1 = nonstationary unit root process\n")
% vratiotest 0 = random walk | 1 = not random walk
fprintf("vratiotest 0 = random walk       | 1 = not random walk\n")
fprintf("--------------------------------\n");
fprintf("Daily\n");
fprintf("----------------\n");
fprintf("Log Price\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(d));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(d));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(d));
fprintf("Log Returns\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(dlr));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(dlr));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(dlr));

fprintf("\nWeekly\n");
fprintf("----------------\n");
fprintf("Log Price\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(w));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(w));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(w));
fprintf("Log Returns\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(wlr));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(wlr));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(wlr));

fprintf("\nMonthly\n");
fprintf("----------------\n");
fprintf("Log Price\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(m));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(m));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(m));
fprintf("Log Returns\n");
fprintf("----------------\n");
fprintf("adftest for      %s: %f\n", name, adftest(mlr));
fprintf("kpsstest for     %s: %f\n", name, kpsstest(mlr));
fprintf("vratiotest for   %s: %f\n", name, vratiotest(mlr));

end

