function [ PPG_START,PK_AMP,PK_IND] = ecg_pkdetect(ecg_sig,Fs)
% ecg_pkdetect computes the indices and amplitude of R-peaks in ECG

%% INPUT
% ecg_sig - input signal
% Fs - Sampling Frequency

%% OUTPUT
% PPG_START  - start sequence for PPG
% PK_AMP - amplitude of R-peak
% PK_IND - Indices of R-peak

%% STAGE -1

% 1.1 BPF Chebyshev - I
% Fs = 125;  % Sampling Frequency

N1      = 6;   % Order
F2_pass1 = 6;   % First Passband Frequency
F2_pass2 = 18;  % Second Passband Frequency
A2_pass  = 1;   % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h1  = fdesign.bandpass('N,Fp1,Fp2,Ap', N1, F2_pass1, F2_pass2, A2_pass, Fs);
Hd1 = design(h1, 'cheby1');

sos1 = Hd1.sosMatrix;
G1 = Hd1.ScaleValues;

ecg_bpf1=filtfilt(sos1,G1,ecg_sig);
% subplot(3,1,2);
% plot(t,ecg_bpf1);


% 1.2 First Order Forward Differencing
a2 =diff(ecg_bpf1);
a3=ecg_bpf1(1)-ecg_bpf1(length(ecg_bpf1));
ecg_diff = [a2;a3];

% 1.3 Amplitude Normalisation
ecg_diff=ecg_diff/max(abs(ecg_diff));

%% STAGE - 2

% 2.1 MEASURES
abs1 = abs(ecg_diff);       % absolute value 
ener1 = ecg_diff.^2;      % Energy value
shan_entrop = -abs1.*log(abs1);     % Shannon entropy
shan_ener =ener1.*log(ener1);     % Shannon Energy

seg_len = 15;
env1= [shan_ener;zeros(seg_len,1)];

for l=1:length(shan_ener)
    env2 = env1(l:(l+seg_len-1));
    env3(l) = (-1/seg_len)*sum(env2);
    mean1=mean(env3);
    sd1=std(env3);
    shan_envel(l)= (env3(l)-mean1);
end

% 2.2 Zero Phase Filtering
% Fs = 125;  % Sampling Frequency
N2  = 30;  % Order
Fc =2;   % Cutoff Frequency
[z1,p1] = butter(N2,0.6,'low');
ecg_zer1=filtfilt(z1,p1,shan_envel');

%% STAGE -3
% 3.1 Hilbert Transformation
h1=hilbert(ecg_zer1);          % Analytic signal
R1=imag(h1);                  % Hilbert Transform of R

% 3.2 Moving Average Filter
N2=312;
q1=0;
m=length(R1);
ma=[R1;zeros(N2,1)];
for i=1:m
    sum1=0;
    k=i;
    for j=k:N2+k
        sum1=sum1+ma(j);
    end 
    
    q1(i)=sum1/N2;
end
q1=q1';
% Hilbert - MA
HMA1=R1-q1;

%% STAGE - 4 
% 4.1 Positive Zero Crossing Detection
r1 = diff(sign(HMA1));
zer_up=find(r1>0);
pk=HMA1(zer_up);

%% STAGE -5

% 5.1 Exact Maximal Sope Detection
Win_l = 25;   PK_AMP = 0; ind = 0;  T_ECG = 0;
pk_1 = [zeros(Win_l,1);ecg_sig;zeros(Win_l,1)];
for i=1:length(zer_up)
   
    rang = zer_up(i):(zer_up(i)+50);
    pk_2 = pk_1(rang);
    [x4(i),I1(i)] = max(pk_2);
    I(i) = rang(I1(i));
    ind (i) = I(i)-Win_l;
    PK_AMP(i) = ecg_sig(I(i)-Win_l);
    PK_IND(i) = ind(i);
%     T_ECG(i) = t(PK_IND(i));
 
   
end
PPG_START = PK_IND(1);

end

