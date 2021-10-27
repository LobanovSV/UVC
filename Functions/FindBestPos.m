function [sco,posS,posE] = FindBestPos(seq,qual,ML,pos0)

Ns = length(seq);

if Ns ~= size(ML,2)
    ML = ML(:,1:Ns);
end

% SaveAlgCoins(seq,ML,'Read')
M = (ML ~= seq) * qual.';

[sco,row] = mink(M,2);

% Consider multi minimums
if sco(1) == sco(2)
    row = find(M == sco(1));
    row = row(randperm(length(row),1));
else
    row = row(1);
end
sco = sco(1);

posS = pos0 - row + 1;
posE = pos0 - row + Ns;