function [ppg_pw] = pts_ppgwidth(morph,P)
% Computing pulse width at different level (5%, 25%, 50%, 75%)

footin = morph.pts_base;
PW05 = []; PW25 = []; PW50 = []; PW75 = [];

for k = 1:length(footin)-1
    PS = P(footin(k):footin(k+1));
        pw50 = pulsewidth(PS);   
        if isempty(pw50)
            pw50 = 0;
        else  pw50 = max(pw50);
        end
        pw05 = pulsewidth(PS,'MidPercentReferenceLevel',5); 
        if isempty(pw05)
            pw05 = 0;
        else  pw05 = max(pw05);
        end

        pw25 = pulsewidth(PS,'MidPercentReferenceLevel',25);
        if isempty(pw25)
            pw25 = 0;
        else pw25 = max(pw25);
        end

        pw75 = pulsewidth(PS,'MidPercentReferenceLevel',75);
        if isempty(pw75)
            pw75 = 0;
        else  pw75 = max(pw75);
        end
        [PW05] = [PW05;pw05]; [PW25] = [PW25;pw25]; [PW50] = [PW50;pw50];
        [PW75] = [PW75;pw75];
end
        ppg_pw = [PW05(:),PW25(:),PW50(:),PW75(:)];

end

