function [Err,Alg,Allele] = CreateFinalAlignment(Allele,SEQ,QUAL, ...
    POSS,PAIR,LEFT,Ind,posW,maxSL,dPos,dPOS,FracReads,Qmax,TextType)

%% Parameters

TolDist = 4;

%% Find errors of final alignment

% Find STRs and expand Alleles
[Allele,posW,STR] = FindAddSTR(Allele,posW);

[PCR,isN,IndN] = FindPCRsREF(Allele,STR,Qmax);

Np = size(Ind,1);
[dEp,iPCRp,POSp] = FindFinalErrors(PCR,STR,SEQ,QUAL,POSS,PAIR,LEFT, ...
    Np,maxSL,dPos,dPOS);

%% Remove alleles with insufficient evidence

dN1 = size(unique(POSp(dEp > 0,[1,4],1),'rows'),1);
dN2 = size(unique(POSp(dEp < 0,[1,4],2),'rows'),1);

% Both alleles are bad
if dN1 == 0 && dN2 == 0
    if sum(Allele{1}==Allele{end}) < sum(Allele{2}==Allele{end})
        dN2 = 1;
    else
        dN1 = 1;
    end
end

if min(dN1,dN2) < FracReads * max(dN1,dN2)
    if dN1 < dN2
        Allele{1} = Allele{2};
        isN{1} = isN{2};
        IndN{1} = IndN{2};
        iPCRp(:,1) = iPCRp(:,2);
        POSp(:,:,1) = POSp(:,:,2);
    end
    Allele{2} = '';
    dEp(:) = inf;
end

%% Sort reads to three Alleles

if dN1 < dN2
    J = [1;2;2;2];
else
    J = [1;2;1;1];
end
JJ = {find(dEp > 0); find(dEp < 0); find(dEp == 0); []};

nSNP = NaN(1,2);
N2 = length(Allele{1});
nE_Allele_p = zeros(1,N2);
for i = 1:3
    I = J(i);
    if isempty(JJ{i})
        continue
    end
    
    % Sort
    [~,jj] = sortrows(POSp(JJ{i},[1,3,2,4],I));
    JJ{i} = JJ{i}(jj);
    
    % Calculate number of SNPs
    if i < 3
        iiAllele = find(Allele{i} ~= Allele{end-1});
        nSNP(i) = ~isempty(iiAllele) + sum(diff(iiAllele)~=1);
        nE_Allele_p = nE_Allele_p + (Allele{i} ~= Allele{end});
    end
end
if isempty(JJ{end-1})
    Allele{end-1} = '';
else
    Allele{end-1} = Allele{I};
end

%% Window positions

ppW = 1:N2;
ppW = ppW >= posW(1) & ppW <= posW(2);

%% Loop

NTT = length(TextType);
Alg = cell(4,1);
nE_Alg_p = zeros(1,N2);
NE_max = zeros(NTT,N2);
sE = zeros(1,N2);
for i = 1:4
    if isempty(Allele{i})
        continue
    end
    I = J(i);
    N1 = length(JJ{i});
    Algt = repmat(blanks(N2),2*N1,1);
    l = 0;
    for j = 1:N1
        ij = JJ{i}(j);
        ipcr = iPCRp(ij,I);
        for k = 1:2
            ijk = Ind(ij,k);
            if ~ijk
                continue
            end
            poss = IndN{I}{ipcr}(POSp(ij,2*k-1,I));
            pose = IndN{I}{ipcr}(POSp(ij,2*k,I));
            
            Nseq = pose + 1 - poss;
            if length(SEQ{ijk}) == Nseq
                seq = SEQ{ijk};
                qual = QUAL{ijk};
            else
                isn = isN{I}(ipcr,poss:pose);
                seq = repmat('*',1,Nseq);
                seq(isn) = SEQ{ijk};
                qual = repmat(Qmax,1,Nseq);
                qual(isn) = QUAL{ijk};
            end
            
            % Remove quality outside of the region of interest
            pp = poss:pose;
            qual(pp<posW(1)|pp>posW(2)) = 0;
            
            % Find line
            l = l + (k == 1 || ~Ind(ij,1) || ...
                poss < IndN{I}{ipcr}(POSp(ij,2,I)) + TolDist);
            
            % Save alignment
            Algt(l,pp) = seq;
            if k == 1
                if poss ~= 1
                    Algt(l,poss-1) = '<';
                end
            else
                if pose ~= N2
                    Algt(l,pose+1) = '>';
                end
            end
            % Calculate error
            jj = seq~=Allele{i}(pp);
            if any(jj)
                ppii = pp(jj);
                qualii = qual(jj);
                nE_Alg_p(ppii) = nE_Alg_p(ppii) + qualii;
                NE_max(:,ppii) = NE_max(:,ppii) + qualii .* ...
                    (seq(jj)==TextType);
            end
            sE(pp) = sE(pp) + qual;
        end
    end
    
    Alg{i} = Algt(1:l,ppW);
    Allele{i} = Allele{i}(ppW);
end

%% Normalize Err and find mean

nE_mean = (1+isempty(Allele{2})) * mean(nE_Allele_p(ppW)) + ...
    mean(nE_Alg_p(ppW)) / Qmax;

NE_max = max(NE_max(:,ppW),[],1);
nE_max = max(NE_max ./ sE(ppW));
NE_max = max(NE_max) / Qmax;

sEwrong = sum(nE_Alg_p(ppW));
sE = sum(sE(ppW));

nSNP = max(nSNP) / sum(Allele{end} ~= '*');

%% Combine errors into the structure

Err = struct('nE_mean',nE_mean,'NE_max',NE_max,'nE_max',nE_max, ...
    'sEwrong',sEwrong,'sE',sE,'nSNP',nSNP,'Qmax',Qmax);

%% Truncate

ppRem = true(1,length(Allele{1}));
for i = 1:4
    if isempty(Allele{i})
        continue
    end
    ppRem = ppRem & (Allele{i}=='*');
    
    AlgRem = Alg{i}==' ' | Alg{i}=='<' | Alg{i}=='>';
    ii = all(AlgRem,2);
    if any(ii)
        Alg{i}(ii,:) = [];
        AlgRem(ii,:) = [];
    end
    ppRem = ppRem & all(AlgRem | Alg{i}=='*',1);
end
if any(ppRem)
    for i = 1:4
        if isempty(Allele{i})
            continue
        end
        Allele{i}(ppRem) = [];
        Alg{i}(:,ppRem) = [];
    end
end