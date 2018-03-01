
function Y = genrateGamCurrent(TotalSimTime)

dt =  0.05;
if nargin < 1;
    TotalSimTime = 100; %<----300
end;

times = 0:dt:TotalSimTime;
Alpha = 3 + 1*rand(1,1); % <----  3, 5
Beta =  10 + 3*rand(1,1); % <---- 15, 18

Y = gampdf(times, Alpha, Beta);
Y = Y./max(Y);