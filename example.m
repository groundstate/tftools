%
% Examples to demonstrate the use of tftools
%

% Generate some 'white' noise
x = rand(1,10000); 
rate = 1.0; % sampling rate in Hz

tau = zeros(1,13);
for i=1:13
    tau(i)=2^(i-1);
end;

% Frequency noise
phase=0;
gaps=0;
% Calculate the overlapping Allan deviation, no gaps
[wfn_a_dev, wfn_a_err, wfn_n_a, wfn_a_new_tau]  = adev(x,rate,tau,1,phase,gaps);

[wfn_m_dev, wfn_m_err, wfn_n_m, wfn_m_new_tau]  = mdev(x,rate,tau,phase);

% Longer averaging times are possible with TOTDEV
% but it's not recommended for t>T/2, where T is the sequence length
tot_tau = [tau 5000 6000 7000 8000 9000];
[wfn_tot_dev, wfn_tot_err, wfn_n_tot, wfn_tot_new_tau] = totdev(x,rate,tot_tau,phase,gaps);

% Phase noise
phase=1;
[wpn_a_dev, wpn_a_err, wpn_a_n, wpn_a_new_tau]  = adev(x,rate,tau,1,phase,gaps);
[wpn_m_dev, wpn_m_err, wpn_n_m, wpn_m_new_tau]  = mdev(x,rate,tau,phase);
[wpn_tot_dev, wpn_tot_err, wpn_n_tot, wpn_tot_new_tau] = totdev(x,rate,tot_tau,phase,gaps);

% Compare
figure(1);
loglog(wfn_a_new_tau,wfn_a_dev,'o-r');
hold on;
loglog(wfn_tot_new_tau,wfn_tot_dev,'+-b');
loglog(wfn_m_new_tau,wfn_m_dev,'+-','Color',[.5 .5 0]);

loglog(wpn_a_new_tau,wpn_a_dev,'*-g');
loglog(wpn_tot_new_tau,wpn_tot_dev,'x-m');
loglog(wpn_m_new_tau,wpn_m_dev,'+-','Color',[0 .7 .7]);

legend('WFN:ADEV','WFN:TOTDEV','WFN:MDEV','WPN:ADEV','WPN:TOTDEV','WPN:MDEV');
xlabel('\tau (s)');
ylabel('deviation');
title('Comparison of ADEV, MDEV and TOTDEV: white phase noise and white frequency noise ');
hold off;

% Now remove some data
tgaps=1:10000; % measurement times
tgaps(1237:1237)=[]; % remove some measurements
xgaps=x;
xgaps(1237:1237)=[];

% Detect and tag gaps
[tgaps_fixed,xgaps_fixed] = markgaps(tgaps,xgaps,1.0);

% display(['Length of data was ' num2str(length(xgaps)) ', now ' num2str(length(xgaps_fixed))]);

gaps=1;
phase=1;
% Compute Allan deviation with tagged data
[gaps_a_dev, gaps_a_err, gaps_a_n, gaps_a_new_tau]  = adev(xgaps_fixed,rate,tau,1,phase,gaps);

% Compare Allan deviation gaps/no gaps
figure(2);
semilogx(wpn_a_new_tau,1.0 - wpn_a_dev ./ gaps_a_dev,'o-');
xlabel('\tau (s)');
ylabel('1 - ADEV(no gap)/ADEV(gap)');
title('Comparison of ADEV with and without gaps - white phase noise');

% Example of use of Theo1
% Generate frequency data with white noise as per NIST algorithm
% (see validation.m)
x = 1:1000;
x(1) = 1234567890;
for i=1:999
    x(i+1)= mod(16807*x(i),2147483647);
end;
x = x / 2147483647;

rate=1.0;
tau = [1 2 4 8 16 32 64 128 256 499];
phase=0;
gaps=0;

[wfn_a_dev, wfn_a_err, wfn_n_a, wfn_a_new_tau]  = adev(x,rate,tau,1,phase,gaps);

% tau = 0.75*(m-1)*tau and m must be even
mTheo1=[ 10 16 32 64 128 256 512 768 998 1000];
%tauTheo1= 0.75*(mTheo1-1); % tau = 1.0
tauTheo1= 0.75*mTheo1;
[wfn_theo1_dev, wfn_theo1_err, wfn_n_theo1, wfn_theo1_new_tau]  = theo1dev(x,rate,tauTheo1,phase,1);
[lowerlim,upperlim]= theo1devci(length(x)+1,mTheo1,68,'white fm'); % add 1 because it's for phase

f3=figure(3);
errorbar(wfn_theo1_new_tau,wfn_theo1_dev,lowerlim .* wfn_theo1_dev, upperlim .* wfn_theo1_dev,'+-');
ax = get(f3,'CurrentAxes');
set(ax,'XScale','log','YScale','log');
hold on;
loglog( wfn_a_new_tau,wfn_a_dev,'O-');
loglog([1.0 1000.0],[1.0/sqrt(12.0) 1.0/sqrt(12)/sqrt(1000)],'-'); 
legend('THEO1DEV','ADEV','expected');
xlabel('\tau (s)');
ylabel('DEV');
title('Comparison of ADEV and THEO1DEV - white frequency noise');

hold off;

%
% Example of use of confidence interval estimates
% Use ADEV calculation from previous example
%
% Limits for 95% confidence
% Here, since the rate equals 1.0, tau=m
[lowerlim,upperlim]= adevci(length(x),tau,95,'white fm');

f4=figure(4);
errorbar(wfn_a_new_tau,wfn_a_dev,lowerlim .* wfn_a_dev, upperlim .* wfn_a_dev,'+-');
ax = get(f4,'CurrentAxes');
set(ax,'XScale','log','YScale','log');
hold on;
loglog([1.0 1000.0],[1.0/sqrt(12.0) 1.0/sqrt(12)/sqrt(1000)],'-'); 
hold off;
xlabel('\tau (s)');
ylabel('ADEV');
title('ADEV with confidence limits');
