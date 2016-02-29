%
% Examples to demonstrate the use of tftools
%
clear all;

% Generate white phase noise
x = rand(1,10000);
rate = 1.0; % sampling rate in Hz

tau = zeros(1,12);
for i=1:12
    tau(i)=2^i;
end;

% Calculate the overlapping Allan deviation, no gaps
[adv, adverr, nadv, newtau]  = adev(x,rate,tau,1,1,0);

% Longer averaging times are possible with TOTDEV
tottau = [tau 5000 6000 7000 8000 9000];
[tdv, tdverr, ntdv, tnewtau] = totdev(x,rate,tottau,1,0);

% Compare
figure(1);
loglog(newtau,adv,'o-');
hold on;
loglog(tnewtau,tdv,'+-');
legend('ADEV','TOTDEV');
xlabel('\tau (s)');
ylabel('deviation');
title('Comparison of ADEV and TOTDEV - white phase noise');
hold off;

% Now remove some data
tgaps=1:10000; % measurement times
tgaps(1237:2679)=[]; % remove some measurements
xgaps=x;
xgaps(1237:2679)=[];

% Detect and tag gaps
[tgapsfixed,xgapsfixed] = markgaps(tgaps,xgaps,1.0);

display(['Length of data was ' num2str(length(xgaps)) ', now ' num2str(length(xgapsfixed))]);

% Compute Allan deviation with tagged data
[gadv, gadverr, gnadv, gnewtau]  = adev(xgapsfixed,rate,tau,1,1,1);

% Compare Allan deviation 
figure(2);
semilogx(newtau,1.0 - adv ./ gadv,'o-');
xlabel('\tau (s)');
ylabel('1 - ADEV(no gap)/ADEV(gap)');
title('Comparison of ADEV with and without gaps - white phase noise');


