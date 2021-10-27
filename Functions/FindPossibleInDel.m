function [INDEL,IIREM] = FindPossibleInDel(REF,SEQ,QUAL,POSS,POSE, ...
    PAIR,LEFT,maxSL,dPos,dPOS,MinReads,Qmax)

%% Parameters

MaxNINDEL = 10;
FracReads = 0.05:0.025:0.2;
MaxSco = 0.1 * maxSL;

if isempty(POSE)
    POSE = POSS;
end

%% InDels

N = length(POSS);
Npair = max(PAIR);
INDEL = cell(N,1);
SCO = zeros(N,1);
Ind = cell(N,1);
PairAlign = false(Npair,1);
[~,iisort] = sort(POSS);
Pos = [0, 0];
for i = 1:N
    j = iisort(i);
    if strcmp(SEQ{j},REF(POSS(j):POSE(j)))
        continue
    end
    PairAlign(PAIR(j)) = true;
    if POSS(j) > Pos(2)
        Pos = [POSS(j), POSS(j)+dPos];
        [M0,pos0] = MatchMatrix(REF,Pos,dPOS,maxSL);
    end
    [~,~,SCO(j),INDEL{j}] = AlignRead(SEQ{j},M0,pos0);
    Ind{j} = repmat(i,size(INDEL{j},1),1);
end

%% Remove reads with bad score

% Both mate pairs must have bad score
IIPAIRREM = true(Npair,1);
ii = PAIR(LEFT);
IIPAIRREM(ii) = IIPAIRREM(ii) & SCO(LEFT) > MaxSco;
ii = PAIR(~LEFT);
IIPAIRREM(ii) = IIPAIRREM(ii) & SCO(~LEFT) > MaxSco;

if any(IIPAIRREM)
    IIREM = ismember(PAIR,find(IIPAIRREM));
    SEQ(IIREM) = [];
    QUAL(IIREM) = [];
    POSS(IIREM) = [];
    PAIR(IIREM) = [];
    LEFT(IIREM) = [];
    INDEL(IIREM) = [];
    Ind(IIREM) = [];
    PairAlign(IIPAIRREM) = [];
    [~,~,PAIR] = unique(PAIR);
else
    IIREM = [];
end

%% Merge and leave InDels with at least MinReads

% Merge
[INDEL,~,iINDEL] = unique(cell2mat(INDEL),'rows');
NINDEL = size(INDEL,1);
if NINDEL == 0
    INDEL = zeros(0,2);
    return
end
[uInd,~,iInd] = unique(cell2mat(Ind));
InDel = false(NINDEL,length(uInd));
InDel(iINDEL+NINDEL*(iInd-1)) = true;
InDel(:,sum(InDel,1)==1) = [];
InDel = unique(InDel.','rows').';
InDel = [false(NINDEL,1),eye(NINDEL,NINDEL,'logical'),InDel];
NInDel = size(InDel,2);

% Remove the same InDels
REFi = cell(NInDel,1);
for i = 1:NInDel
    REFi{i} = AddInDel(REF,INDEL(InDel(:,i),:));
end
[~,ii] = unique(REFi,'stable');
if length(ii) < NInDel
    InDel = InDel(:,ii);
    ii = ~any(InDel,2);
    if any(ii)
        INDEL(ii,:) = [];
        if isempty(INDEL)
            return
        end
        InDel(ii,:) = [];
    end
end

% Skip pairs with two perfectly matched reads
PairAlign = PairAlign(PAIR);
SEQ = SEQ(PairAlign);
QUAL = QUAL(PairAlign);
POSS = POSS(PairAlign);
[~,~,PAIR] = unique(PAIR(PairAlign));
LEFT = LEFT(PairAlign);

% Find errors and positions of pairs
for i = 1:length(FracReads)
    [~,~,iRem] = SimplePairAlignment(REF,SEQ,QUAL,POSS,PAIR,LEFT, ...
        INDEL,InDel,maxSL,dPos,dPOS,MinReads,Qmax,FracReads(i));
    
    InDel(:,iRem) = [];
    ii = ~any(InDel,2);
    if any(ii)
        INDEL(ii,:) = [];
        InDel(ii,:) = [];
    end
    
    if size(INDEL,1) < MaxNINDEL
        break
    end
end