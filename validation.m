
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
