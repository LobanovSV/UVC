function [El,posSl,posEl,Er,posSr,posEr] = FindBestPairPos( ...
    pos0l,seql,quall,pos0r,seqr,qualr,REF,dPOS,MinSL)

Nsl = length(seql);
Nsr = length(seqr);

[ML,pos0l] = MatchMatrix(REF,pos0l,dPOS,Nsl);
[MR,pos0r] = MatchMatrix(REF,pos0r,dPOS,Nsr);

Ml = (ML ~= seql) * quall.';
Mr = (MR ~= seqr) * qualr.';

% Check
if dPOS < 10
    error('Change dPOS')
end

[El,rowl] = mink(Ml,10);
[Er,rowr] = mink(Mr,10);

EE = El + Er.';
[~,ii] = sort(EE(:));
for i = 1:length(ii)
    [il,ir] = ind2sub([10,10],ii(i));
    posSl = pos0l - rowl(il) + 1;
    posEr = pos0r - rowr(ir) + Nsr;
    if posEr >= posSl + MinSL
        posEl = pos0l - rowl(il) + Nsl;
        posSr = pos0r - rowr(ir) + 1;
        El = El(il);
        Er = Er(ir);
        return
    end
end

%% Do not found

[il,ir] = ind2sub([10,10],ii(1));
posSl = pos0l - rowl(il) + 1;
posEr = pos0r - rowr(ir) + Nsr;
posEl = pos0l - rowl(il) + Nsl;
posSr = pos0r - rowr(ir) + 1;
El = El(il);
Er = Er(ir);