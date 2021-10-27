function [Ep,POSp,iRem,isPCRp] = SimplePairAlignment(REF,SEQ,QUAL, ...
    POSS,PAIR,LEFT,INDEL,InDel,maxSL,dPos,dPOS,MinReads,Qmax,FracReads)
% Find errors and positions of pairs. Additionally find INDELs which should
% be removed

%% Input parameters

N1 = size(INDEL,1);
OnlyINDEL = isempty(InDel);
if OnlyINDEL
    InDel = [false(N1,1),eye(N1,N1,'logical')];
elseif any(InDel(:,1))
    error('First column must be false')
end

%% REF errors

Np = max(PAIR);
N2 = size(InDel,2);
Ep = zeros(Np,N2);
POSp = repmat(inf(Np,1) .* [-1,-1,1,1],1,1,N2);
isPCRp = false(Np,N2);
RIGHT = ~LEFT;

[E,POSS,POSE,isPCR] = SimpleAlignmentSTR(REF,SEQ,QUAL,POSS,maxSL, ...
    dPos,dPOS,Qmax,[],PAIR,LEFT);

ii = PAIR(LEFT);
Ep(ii,1) = E(LEFT);
isPCRp(ii,1) = isPCR(LEFT);
POSp(ii,1,1) = POSS(LEFT);
POSp(ii,2,1) = POSE(LEFT);
ii = PAIR(RIGHT);
Ep(ii,1) = Ep(ii,1) + E(RIGHT);
isPCRp(ii,1) = isPCRp(ii,1) | isPCR(RIGHT);
POSp(ii,3,1) = POSS(RIGHT);
POSp(ii,4,1) = POSE(RIGHT);
if isempty(INDEL)
    iRem = false;
    return
end

%% Find INDEL Shift and ReCalc

[pS,pE] = INDEL2POS(INDEL);
pS = pS.';
pE = pE.';
dp = 1 - diff(INDEL,1,2).';
Shift = (POSS > pE) .* dp;
ReCalc = (POSS <= pE + dPOS) & (POSE >= pS - dPOS);

%% ALT errors

for i = 2:N2
    jj = InDel(:,i);
    REFi = AddInDel(REF,INDEL(jj,:));
    shift = sum(Shift(:,jj),2);
    Ei = E;
    isPCRi = isPCR;
    POSSi = POSS + shift;
    POSEi = POSE + shift;
    ii = any(ReCalc(:,jj),2);
    [Ei(ii),POSSi(ii),POSEi(ii),isPCRi(ii)] = ...
        SimpleAlignmentSTR(REFi,SEQ(ii),QUAL(ii),POSSi(ii),maxSL,dPos, ...
        dPOS,Qmax,[],PAIR(ii),LEFT(ii));
    ii = PAIR(LEFT);
    Ep(ii,i) = Ei(LEFT);
    isPCRp(ii,i) = isPCRi(LEFT);
    POSp(ii,1,i) = POSSi(LEFT);
    POSp(ii,2,i) = POSEi(LEFT);
    ii = PAIR(RIGHT);
    Ep(ii,i) = Ep(ii,i) + Ei(RIGHT);
    isPCRp(ii,i) = isPCRp(ii,i) | isPCRi(RIGHT);
    POSp(ii,3,i) = POSSi(RIGHT);
    POSp(ii,4,i) = POSEi(RIGHT);
end
if nargout < 3
    return
end

%% Find INDELs which should be removed

% Calculate INDEL intersection matrix
InterSect = pE.' >= pS & pS.' <= pE;

iRem = false(1,N2);
for i = 2:N2
    if sum(~isPCRp(:,i) & Ep(:,i) < Ep(:,1)) < MinReads
        continue
    end
    Ni = size(unique(POSp(~isPCRp(:,i) & Ep(:,i) < Ep(:,1), ...
        [1,4],i),'rows'),1);
    N1 = size(unique(POSp(~isPCRp(:,1) & Ep(:,1) < Ep(:,i), ...
        [1,4],1),'rows'),1);
    iRem(i) = Ni >= max(MinReads,FracReads*N1);
    if iRem(i)
        for j = find(iRem)
            if i == j
                continue
            end
            Ni = size(unique(POSp(~isPCRp(:,i) & Ep(:,i) < Ep(:,j), ...
                [1,4],i),'rows'),1);
            Nj = size(unique(POSp(~isPCRp(:,j) & Ep(:,j) < Ep(:,i), ...
                [1,4],j),'rows'),1);
            IS = any(InterSect(InDel(:,i),InDel(:,j)),'all');
            if (Ni < MinReads && Nj > Ni + MinReads) || ...
                    (IS && Ni < FracReads*Nj)
                iRem(i) = false;
                break
            elseif (Nj < MinReads && Ni > Nj + MinReads) || ...
                    (IS && Nj < FracReads*Ni)
                iRem(j) = false;
                continue
            end
        end
    end
end

iRem = ~iRem;
if OnlyINDEL
    iRem(1) = [];
else
    iRem(1) = false;
end