function [Allele,Alg,Err,IIREM] = SingleSampleAlignment(REF,SEQ,QUAL, ...
    POSS,POSE,PAIR,LEFT,posW,dPos,dPOS,RAM,FILE,i,AlleleOnly)


%% Default parameters

if nargin < 11 || isempty(RAM)
    RAM = true;
end

%% Parameters

MinReads = 3;
FracReads = 0.05;
M2ProdTol = 1.1;
M2AddTol = 5;

%% Default output

Allele = [];
Alg = [];
Err = struct('nE_mean',0,'NE_max',0,'nE_max',0,'sEwrong',0,'sE',0, ...
    'nSNP',0,'Qmax',0);

%% Load data

if ~RAM
    File = fullfile(DirPCorUNIX,'BAMs','Temp',FILE,[int2str(i) '.mat']);
    load(File,'SEQ','QUAL','PAIR','LEFT','POSS','POSE')
    save(File,'Allele','Alg','Err','-append')
end

%% No reads

if isempty(SEQ)
    IIREM = [];
    return
end

%% Maximum quality

Qmax = max(cellfun(@max,QUAL));
maxSL = 1 + max(POSE - POSS);

%% Find Possible InDels

[INDEL,IIREM] = FindPossibleInDel(REF,SEQ,QUAL,POSS,POSE,PAIR,LEFT, ...
    maxSL,dPos,dPOS,MinReads,Qmax);

%% Remove reads with bad score

if ~isempty(IIREM)
    SEQ(IIREM) = [];
    QUAL(IIREM) = [];
    POSS(IIREM) = [];
    POSE(IIREM) = [];
    PAIR(IIREM) = [];
    [~,~,PAIR] = unique(PAIR);
    LEFT(IIREM) = [];
    if ~RAM
        save(File,'SEQ','QUAL','PAIR','LEFT','POSS','POSE','Allele', ...
            'Alg','Err')
    end
    % No reads
    if isempty(SEQ)
        IIREM = [];
        return
    end
end

%% Find all possible alleles

InDel = FindAllPossibleAlleles(INDEL,REF);

%% Find errors and positions of possible InDels

[Ep,POSp,iRem,isPCRp] = SimplePairAlignment(REF,SEQ,QUAL,POSS,PAIR, ...
    LEFT,INDEL,InDel,maxSL,dPos,dPOS,MinReads,Qmax,FracReads);
if any(iRem)
    InDel(:,iRem) = [];
    Ep(:,iRem) = [];
    POSp(:,:,iRem) = [];
    isPCRp(:,iRem) = [];
end

%% Mismatch for pairs

N = size(Ep,2);
M2 = NaN(N,N);
for i = 1:N
    for j = 1:i
        dN1 = size(unique(POSp(~isPCRp(:,i) & Ep(:,i) < Ep(:,j), ...
            [1,4],i),'rows'),1);
        dN2 = size(unique(POSp(~isPCRp(:,j) & Ep(:,j) < Ep(:,i), ...
            [1,4],j),'rows'),1);
        if min(dN1,dN2) >= FracReads * (dN1 + dN2)
            M2(i,j) = sum(min(Ep(:,i),Ep(:,j)));
        end
    end
end

%% Find minimum possible 

M2Thr = min(M2,[],'all') * M2ProdTol + M2AddTol * Qmax;

ii = (min(M2,[],1) <= M2Thr) | (min(M2,[],2).' <= M2Thr);
if ~all(ii)
    InDel = InDel(:,ii);
    Ep = Ep(:,ii);
    POSp = POSp(:,:,ii);
    isPCRp = isPCRp(:,ii);
    M2 = M2(ii,ii);
    
    ii = any(InDel,2);
    INDEL = INDEL(ii,:);
    InDel = InDel(ii,:);
end

[~,ii] = sort(M2(:));
ii = ii(M2(ii)<=M2Thr).';
[i1,i2] = ind2sub(size(M2),ii);

%% Find indexes of pairs

Ind = zeros(size(Ep,1),2);
Indt = 1:length(PAIR);
Ind(PAIR(LEFT),1) = Indt(LEFT);
Ind(PAIR(~LEFT),2) = Indt(~LEFT);

%% Find two best alleles

Err.nE_mean = inf;
for i = 1:length(i1)
    ii = [i1(i),i2(i)];
    dEp = diff(Ep(:,ii),1,2);
    [Allelet,Algt,Errt] = FindAlignment(INDEL,InDel(:,ii),...
        REF,SEQ,QUAL,POSS,PAIR,LEFT,dEp,POSp(:,:,ii),isPCRp(:,ii), ...
        Ind,posW,maxSL,dPos,dPOS,Qmax,MinReads,FracReads);
    if Errt.nE_mean < Err.nE_mean
        Err = Errt;
        Allele = Allelet;
        Alg = Algt;
    end
end

%% Save

if ~RAM
    if AlleleOnly
        Alg = [];
        Err = [];
    end
    save(File,'Allele','Alg','Err','-append')
    Allele = [];
    Alg = [];
    Err = [];
end