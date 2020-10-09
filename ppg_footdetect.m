function [footin,footamp,pin,pamp] = ppg_footdetect(ppg,Fs)
%MODIFIED PEAK POINTS &FOOT POINT OF EACH PPG PULSE
%% INPUT
% ppg - input ppg signal
% Fs - Sampling Frequency

%% OUTPUT
% footin - Indices of ppg-foot
% footamp - Indices of ppg-peak
% pamp - Amplitude of ppg-peak
% pin - Indices of ppg-peak

if Fs<=125
    [pin,pamp] = ppg_pkdetect(ppg,Fs);
    thresh = 30;
else 
    [pin,pamp] = ppg_pkdetect1000(ppg);
    thresh = 65;
end

 pin(end)=[]; pamp(end)=[]; 
 
 diff1 = diff(pin);       
mean1 = mean(diff1);  std1 = std(diff1);
I1 = find(diff1 < 0.8*(mean1-std1));
        pin(I1+1)=[];   pamp(I1+1)=[];
x = find(diff(pin)<50);
pin(x+1)=[];
pamp(x+1)=[];

% To Find the PPG Base point with Interval range of Two ECG Peaks
if pin(1)>thresh
    for i = 1:length(pin)
        rng2 = [pin(i)-thresh:pin(i)];
        [Base] = ppg(rng2);
        [max2(i),indmax1(i)] = min(Base);

        ind_abase(i) = rng2(indmax1(i));
        ppg_abase(i) = max2(i);
        ind_base(i)  = ind_abase(i);
        ppg_base(i)  = ppg(ind_abase(i));
    end
else
    for i = 1:length(pin)-1
        rng2 = [pin(i+1)-thresh:pin(i+1)];
        [Base] = ppg(rng2);
        [max2(i),indmax1(i)] = min(Base);

        ind_abase(i) = rng2(indmax1(i));
        ppg_abase(i) = max2(i);
        ind_base(i)  = ind_abase(i);
        ppg_base(i)  = ppg(ind_abase(i));
    end
    pin(1) = []; pamp(1) = []; 
end
    footin = [ind_base];
    footacc = [ppg_abase];
    footamp = [ppg_base];
end

