function [row,colL,colR,sco] = AlignReadLMR0(L0Del,L0Ins,M,Row,ColL, ...
    ColR,LL,k)

[row,colL,colR,sco] = AlignReadLMR(L0Del,L0Ins,M,Row,ColL,ColR,LL,k);

% Leave only insertions or deletions
if any(diff(row) > 0) && any(diff(row) < 0)
    [row,colL,colR,sco] = AlignReadLMR(L0Del,inf,M,Row,ColL,ColR,LL,k);
    [rowt,colLt,colRt,scot] = AlignReadLMR(inf,L0Ins,M,Row,ColL,ColR,LL,k);
    if scot < sco
        row = rowt;
        colL = colLt;
        colR = colRt;
        sco = scot;
    end
end