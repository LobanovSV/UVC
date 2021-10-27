function InDel = FindAllPossibleAlleles(INDEL,REF)

% Find unique deletions
N1 = size(INDEL,1);
if N1 == 0
    InDel = false;
    return
end

%% Find start and end

[pS,pE] = INDEL2POS(INDEL);

%% Find unique positions

InDel = [false,true];
Ndel = [0,1]; % Number of deletions
for i = 2:N1
    Intersect = pS(1:i-1) <= pE(i) & pE(1:i-1) >= pS(i);
    ii = ~any(InDel(Intersect,:),1);
    indel = InDel(:,ii);
    N0 = size(InDel,2);
    N1 = size(indel,2);
    Ndel = [Ndel, Ndel(ii) + 1 - sum(indel & (pS(1:i-1) == pE(i) - 1 | ...
        pE(1:i-1) == pS(i) + 1),1)]; %#ok
    InDel = [InDel,indel;false(1,N0),true(1,N1)]; %#ok
end

% Sort with respect to the number of InDdels
[~,ii] = sort(Ndel);
InDel = InDel(:,ii);

%% Remove alleles with the same sequence

N = size(InDel,2);
pS = min(pS);
pE = max(pE);
REF = REF(pS:pE);
INDEL = INDEL - (pS - 1);
REFi = cell(N,1);
for i = 1:N
    REFi{i} = AddInDel(REF,INDEL(InDel(:,i),:));
end
[~,ii] = unique(REFi,'stable');
InDel = InDel(:,ii);