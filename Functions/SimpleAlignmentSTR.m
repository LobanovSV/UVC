function [E,POSS,POSE,isPCR] = SimpleAlignmentSTR(REF,SEQ,QUAL, ...
    POSS,maxSL,dPos,dPOS,Qmax,pos,PAIR,LEFT)


%% Default parameters

if nargin < 9
    pos = [];
end
if nargin < 10
    PAIR = [];
end
if nargin < 11
    LEFT = [];
end

%% Parameters

%% Simple alignment

[E,POSS,POSE] = SimpleAlignment(REF,SEQ,QUAL,POSS,maxSL,dPos, ...
    dPOS,PAIR,LEFT);

%% Find STRs

isPCR = false(size(POSS));
[~,ML,PosS,PosE] = FindSTR(REF,pos);
if isempty(ML)
    return
end

%% Correct errors of PCR dublicates

for i = 1:length(ML)
    ml = ML(i);
    posS = PosS(i);
    posE = PosE(i);
    NR = (posE + 1 - posS) / ml;
    PosM = posS - 1 + round(NR/2) * ml;
    
    [NPCR,Qmaxi] = FindPCRerr(NR,Qmax,ml);
    
    % Choose intersecting reads
    II = (POSS < posS + ml * NPCR & POSE > posE) | ...
        (POSS < posS & POSE > posE - ml * NPCR);
    if ~any(II)
        continue
    end
    Ei = E(II);
    SEQi = SEQ(II);
    QUALi = QUAL(II);
    POSSi = POSS(II);
    POSEi = POSE(II);
    isPCRi = isPCR(II);
    PAIRi = PAIR(II);
    LEFTi = LEFT(II);
    
    for j = [-NPCR:-1,1:NPCR]
        REFj = REF([1:PosM-1+j*ml,PosM:end]);
        [Et,POSSt,POSEt] = SimpleAlignment(REFj,SEQi,QUALi,POSSi, ...
            maxSL,dPos,dPOS,PAIRi,LEFTi);
        Et = Et + abs(j) * Qmaxi;
        
        ii = Et < Ei;
        if any(ii)
            Ei(ii) = Et(ii);
            isPCRi(ii) = true;
            
            posst = POSSt(ii);
            poset = POSEt(ii);
            jj = posst >= posS;
            posst(jj) = max(posS,posst(jj) - j*ml);
            jj = poset >= posS;
            poset(jj) = poset(jj) - j*ml;
            POSSi(ii) = posst;
            POSEi(ii) = poset;
        end
    end
    
    E(II) = Ei;
    POSS(II) = POSSi;
    POSE(II) = POSEi;
    isPCR(II) = isPCRi;
end