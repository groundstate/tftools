% Test Theo1 deviation bias correction

% Generate some 'white' noise
x = rand(1,10000); 
rate = 1.0; % sampling rate in Hz
phase =1;
gaps =0;
overlapping=1;

tau = zeros(1,13);
for i=1:13
    tau(i)=2^(i-1);
end;

% Allan deviation
[wpn_a_dev, wpn_a_err, wpn_n_a, wpn_a_new_tau]  = adev(x,rate,tau,overlapping,phase,gaps);


% Theo1
mtheo1 = [tau(5:13) 8192 9950];
tau_theo1 = 0.75*mtheo1;
[wpn_theo1_dev, wpn_theo1_err, wpn_n_theo1, wpn_theo1_new_tau]  = theo1dev(x,rate,tau_theo1,phase);

% Do the bias correction
[ bias, wpn_theo1_dev_corr ] = theo1devbias(mtheo1,'white pm',wpn_theo1_dev);

loglog(wpn_a_new_tau,wpn_a_dev,'+-');
hold on;
loglog(wpn_theo1_new_tau,wpn_theo1_dev,'o-');
loglog(wpn_theo1_new_tau,wpn_theo1_dev_corr,'*-');
legend('ADEV','THEO1','THEO1 (corrected)');
xlabel('\tau (s)');
ylabel('DEV');
title('Correction of THE01 deviation - white PM');
hold off;