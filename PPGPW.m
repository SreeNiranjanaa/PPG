function [pw] = PPGPW(morph,ppg)
%   Manual Calculation of Pulse width
%   Input - ppg (resampled PPG)
%           morph - points at various parts of PPG
%  Output - pw50 (Pulse width at 50%)


%%
% Butterworth filter design
[b,a] = butter(4,0.009,'low');
% Filtering using zero-phase filtering
ppg_sig = filter(b,a,ppg);

footin = morph.pts_base;
bs1 = footin(1:end-1);
bs2 = footin(2:end);
for k=1:length(bs1)
    seg1 = bs1(k):bs2(k);
    seg2 = ppg_sig(seg1);
    seg2 = (seg2-min(seg2))/(max(seg2)-min(seg2));
    [ipk,~] = find(seg2==1);
    rng1 = 1:ipk;
    rng2 = ipk:length(seg2);
    seg3 = seg2(rng1);
    seg4 = seg2(rng2);
    
    [ibs125,~] = find(seg3 >=0.24 & seg3 <=0.25);
    [ibs225,~] = find(seg4 >=0.24 & seg4 <=0.26);
    ibs225 = rng2(ibs225);
    if (numel(ibs125) ==0 || numel(ibs225) == 0)
        pw25(k) = 0;
    else
        pw25(k) = ibs225(1)-ibs125(end);
    end
    
    [ibs150,~] = find(seg3 >=0.49 & seg3 <=0.50);
    [ibs250,~] = find(seg4 >=0.49 & seg4 <=0.50);
    ibs250 = rng2(ibs250);
    if (numel(ibs150)==0 || numel(ibs250) == 0)
        pw50(k) = 0;
    else
        pw50(k) = ibs250(end)-ibs150(end);
    end

    [ibs175,~] = find(seg3 >=0.74 & seg3 <=0.75);
    [ibs275,~] = find(seg4 >=0.74 & seg4 <=0.75);
    ibs275 = rng2(ibs275);
    if (numel(ibs175) ==0 || numel(ibs275) == 0)
        pw75(k) = 0;
    else
        pw75(k) = ibs275(1)-ibs175(end);
    end
    
end

pw25 = pw25(:);  
pw50 = pw50(:);  pw75=pw75(:);
pw = [pw25,pw50,pw75];

end

