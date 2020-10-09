function [ SBP,DBP,PK_IND,DBPIND ] = abp_pkdetect(abp_sig,Fs)
% abp_pkdetect computes the indices and amplitude of Systolic
% Blood pressure and Diastolic Blood pressure 

%% INPUT
% abp_sig - input abp signal

%% OUTPUT
% SBP - Systolic Blood pressure
% DBP - Diastolic Blood Pressure
% PK_IND - Indices of SBP
% DBPIND - Indices of DBP

%%
% Fs = 125;  % Sampling Frequency

N1      = 4;   % Order
F2_pass1 = 0.2;   % First Passband Frequency
F2_pass2 = 4.7;  % Second Passband Frequency
A2_pass  = 1;   % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h1  = fdesign.bandpass('N,Fp1,Fp2,Ap', N1, F2_pass1, F2_pass2, A2_pass, Fs);
Hd2 = design(h1, 'cheby1');

sos1 = Hd2.sosMatrix;
G1 = Hd2.ScaleValues;

ppg_f1=filtfilt(sos1,G1,abp_sig);

% acceleration PPG
appg = derivative(ppg_f1,2);

%% Stage -2
% 2.1 MEASURES
abs1 = abs(appg);       % absolute value 
ener1 = appg.^2;      % Energy value
shan_entrop = -abs1.*log(abs1);     % Shannon entropy
shan_ener =-ener1.*log(ener1);     % Shannon Energy

seg_len = 15;
env1= [shan_ener;zeros(seg_len,1)];

for l=1:length(shan_ener)
    env2 = env1(l:(l+seg_len-1));
    env3(l) = (-1/seg_len)*sum(env2);
    mean1=mean(env3);
    sd1=std(env3);
    shan_envel(l)= (env3(l)-mean1);
end

% subplot(2,1,1);plot(ppg_sig);
% subplot(2,1,2);plot(appg);


% 2.2 Zero Phase Filtering
% Fs = 125;  % Sampling Frequency

n2  = 18 ;  % Order
Fc =8;   % Cutoff Frequency

Wn2 = (Fc*2)/Fs;
[z,p] = butter(n2,Wn2,'low');
ppg_zer1=filtfilt(z,p,shan_envel');

% figure();
% subplot(3,1,1);
% plot(shan_ener);
% subplot(3,1,2);  plot(-shan_envel); 
% subplot(3,1,3); plot(-ppg_zer1);

%% STAGE -3
% 3.1 Hilbert Transformation
h=hilbert(-ppg_zer1);          % Analytic signal
R=imag(h);                  % Hilbert Transform of R

% 3.2 Moving Average Filter
N2=500;
q=0;
m=length(R); sum2 = 0;
ma=[R;zeros(N2,1)];
for i=1:m
    sum2=0;
    k=i;
    for j=k:N2+k
        sum2=sum2+ma(j);
    end 
    
    q(i)=sum2/N2;
end
q=q';
% Hilbert - MA
HMA=R-q;

% figure();
% subplot(3,1,1);plot(R);
% subplot(3,1,2); plot(q);
% subplot(3,1,3); plot(HMA);

%% STAGE - 4 
% 4.1 Positive Zero Crossing Detection
r3 = diff(sign(HMA));
zer_up=find(r3>0);
pk=HMA(zer_up);

% figure();
% plot(HMA);
% hold on;
% plot(zer_up,pk,'k^','markerfacecolor',[0 1 0]);

PK_AMP=abp_sig(zer_up);

% figure();
% plot(abp_sig); hold on;
% plot(zer_up,PK_AMP,'ko','markerfacecolor',[1 0 0]);

%% STAGE -5

% 5.1 Exact Maximal Sope Detection
Win_l = 30;   ppg_pk = 0; ind_p = 0;  t_ecg = 0;
pk_1 = [zeros(Win_l,1);abp_sig;zeros(Win_l,1)];
if zer_up(end)<310
    for i=1:length(zer_up)
        rang = zer_up(i):(zer_up(i)+65);
        pk_2 = pk_1(rang);
        [x4(i),I1(i)] = max(pk_2);
        I(i) = rang(I1(i));
        ind_p (i) = I(i)-Win_l;
        ppg_pk(i) = abp_sig(I(i)-Win_l);
    %     t_pk = t(ind_p);  
    end
else
    for i=1:length(zer_up)-1
        rang = zer_up(i):(zer_up(i)+65);
        pk_2 = pk_1(rang);
        [x4(i),I1(i)] = max(pk_2);
        I(i) = rang(I1(i));
        ind_p (i) = I(i)-Win_l;
        ppg_pk(i) = abp_sig(I(i)-Win_l);
    %     t_pk = t(ind_p);  
    end
end


a=find(diff(ind_p)<50); 
temp1=ind_p; 
temp2=ppg_pk;

temp3 = temp1; temp4 = temp2;

if ~isempty(a)
    for i=1:numel(a)
        arr = [a(i),a(i)+1];
        [pk1(i),pk2(i)]=min(temp2(arr));
        tem(i) = arr(pk2(i));
    end
    temp3(tem) = [];
    temp4(tem)=[];
else
    temp3 = temp1; temp4 = temp2;
end

PK_IND = temp3; 
PK_AMP = temp4;

% SBP = (median(PK_AMP));
SBP = PK_AMP;

% For diastolic BP
if numel(temp3)==1
    [drng] = temp3(1):375;
    [damp] = abp_sig(drng);
    [din1,dind1] = min(damp);
    dind12 = drng(dind1);
else
    for m=1:length(temp3)-1
        [drng] = temp3(m):temp3(m+1);
        [damp] = abp_sig(drng);
        [din1(m),dind1(m)] = min(damp);
        [dind12(m)] = drng(dind1(m));
    end
end
% % For last value
% drng2 = temp3(end):375;
% damp2 = ABPsig(drng2);
% [din2,dind2] = min(damp2);
% dind22 = drng2(dind2);

[DBPIND] = [dind12];
[FT_AMP] = [din1];

% DBP = median(FT_AMP);
DBP = FT_AMP;
end

