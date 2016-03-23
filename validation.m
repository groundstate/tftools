
%

format long;

% Expected value from NIST NIST SP1065 p.108
adevNIST  = [2.922319e-01 9.965736e-02 3.897804e-02];
oadevNIST = [2.922319e-01 9.159953e-02 3.241343e-02];
mdevNIST  = [2.922319e-01 6.172376e-02 2.170921e-02];
totdevNIST= [2.922319e-01 9.134743e-02 3.406530e-02];
tdevNIST  = [1.687202e-01 3.563623e-01 1.253382];

% Generate frequency data with white noise as per NIST algorithm
n = 1:1000;
n(1) = 1234567890;
for i=1:999
    n(i+1)= mod(16807*n(i),2147483647);
end;
n = n / 2147483647;

x=freq2phase(n,1.0);
tau=[1 10 100];
rate=1.0;

[adv,  adverr,  nadv]=adev(x,rate,tau,0);  % non-overlapping ADEV
[oadv, oadverr, noadv]=adev(x,rate,tau); % overlapping ADEV
[mdv,  mdverr,  nmdv]=mdev(x,rate,tau);    % 
[tdv,  tdverr,  ntdv]=totdev(x,rate,tau);
[timedev, timedeverr,ntimedev]=tdev(x,rate,tau);

display('Ratio of computed statistics with NIST test data: NIST SP1065 p.108');
display('-------------------------------------------------------------');
display( '             T=1       T=10     T=100');

r = adevNIST ./ adv;
display(['ADEV    ' num2str(r(1),10) ' ' num2str(r(2),10) ' ' num2str(r(3),10)]);

r = oadevNIST ./ oadv;
display(['OADEV   ' num2str(r(1),10) ' ' num2str(r(2),10) ' ' num2str(r(3),10)]);

r = mdevNIST ./ mdv;
display(['MDEV    ' num2str(r(1),10) ' ' num2str(r(2),10) ' ' num2str(r(3),10)]);

r = totdevNIST ./ tdv;
display(['TOTDEV  ' num2str(r(1),10) ' ' num2str(r(2),10) ' ' num2str(r(3),10)]);

r = tdevNIST ./ timedev;
display(['TDEV    ' num2str(r(1),10) ' ' num2str(r(2),10) ' ' num2str(r(3),10)]);

display('-------------------------------------------------------------');
display('Note: NIST values have 7 significant digits.');

% Validation against examples in IEEE 1139-1999
% These are useful because they are short sets and so expose a different
% set of potential bugs

display(' ');
display('Validation against IEEE 1139-1999');
display('-------------------------------------------------------------');

% ADEV Appendix C3
% Expected value for tau=2 is 3.95E-6
x=[0 43.6 89.7 121.6 163.7 208.4 248 289 319.8]; % in us
tau=[1 2];
rate=1.0;
[adv,  adverr,  nadv]=adev(x,rate,tau);  
display(['OADEV (tau=2) : expect  ' num2str(3.95E-6) ', calculate ' num2str(adv(2)*1.0E-6)]);


% MDEV Appendix C3
% Expected value for tau=2 is 2.47E-6
x=[0 43.6 89.7 121.6 163.7 208.4 248 289 319.8]; % in us
tau=[1 2];
rate=1.0;
[mdv,  mdverr,  nmdv]=mdev(x,rate,tau);  
display(['MDEV  (tau=2) : expect  ' num2str(2.47E-6) ', calculate ' num2str(mdv(2)*1.0E-6)]);

% TOTDEV Appendix C4
% Expected value for tau=2 is 1.79E-9
x=[1.08 0.5 2.2 4.68 3.29]; % all in ns
tau=[1 2];
rate=1.0;
[tdv,  tdverr,  ntdv]=totdev(x,rate,tau);
display(['TOTDEV(tau=2) : expect  ' num2str(1.79E-9) ', calculate ' num2str(tdv(2)*1.0E-9)]);
