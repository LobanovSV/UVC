function [A,AisN,G,AMeshW] = FindAandG(RS,RE,Allele,FPOS,isBad,FASTA, ...
    posW,RAM,FILE)
% Allele = {Allele1,Allele2,PCRs,Allele with highest coverage,FASTA}
% iA - FASTA allele index

%% Find unique alleles

N = length(Allele);
A = cell(2*N+1,1);
A{end} = FASTA(posW(1)-1+RS:posW(1)-3+RE);
for i = 1:N
    [Allelei,FPOSi] = LoadAlleleFPOS(Allele,FPOS,RAM,FILE,i);
    for j = 1:2
        if isBad(i,j)
            A{2*i-2+j} = A{end};
        elseif ~isempty(Allelei{j})
            allelet = Allelei{j}(FPOSi(RS)+1:FPOSi(RE)-1);
            allelet(allelet=='*') = [];
            A{2*i-2+j} = allelet;
        elseif j == 2
            A{2*i} = A{2*i-1};
        else
            error('Impossible')
        end
    end
end
[A,ii,G] = unique(A);

%% Find AisN

NA = length(A);
AisN = zeros(NA,RE-RS);
for k = 1:NA
    if G(end) == k
        AisN(k,1:end-1) = 1;
    else
        i = ii(k);
        j = 1 + mod(i-1,2);
        i = 1 + (i - j) / 2;
        [Allelei,FPOSi] = LoadAlleleFPOS(Allele,FPOS,RAM,FILE,i);
        kk = FPOSi(RS)+1:FPOSi(RE)-1;
        sN1 = cumsum([false,Allelei{j}(kk) ~= '*',false]);
        i2 = find(Allelei{end}(kk) ~= '*') + 1;
        AisN(k,:) = sN1([i2,end])-sN1([1,i2]);
    end
end

%% Reorder

iAt = G(end);
ii = [iAt,1:iAt-1,iAt+1:NA];
A = A(ii);
AisN = AisN(ii,:);
G(G==iAt) = 0;
ii = G < iAt;
G(ii) = G(ii) + 1;

%% Convert

G = reshape(G(1:end-1),2,[]).';

%% Find genotype of bad alleles

for j = 1:2
    for i = find(isBad(:,j)).'
        [Allelei,FPOSi] = LoadAlleleFPOS(Allele,FPOS,RAM,FILE,i);
        allelet = Allelei{j}(FPOSi(RS)+1:FPOSi(RE)-1);
        allelet(allelet=='*') = [];
        k = find(strcmp(A,allelet),1);
        if ~isempty(k)
            G(i,j) = k;
        end
    end
end

%% Find AMeshW

AMeshW = eye(NA);

% Full alleles
AFASTA = cell(NA,1);
pS = posW(1)-2+RS;
for i = 1:NA
    AFASTA{i} = [FASTA(1:pS),A{i},FASTA(pS-RS+RE:end)];
end

% Find STRs
for i = 1:NA
    pE = pS + 1 + length(A{i});
    [~,ML,PosS,PosE] = FindSTR(AFASTA{i});
    ii = PosS > pE | PosE < pS;
    ML(ii) = [];
    PosS(ii) = [];
    PosE(ii) = [];
    for j = 1:length(ML)
        ml = ML(j);
        posS = PosS(j);
        posE = PosE(j);
        NR = (posE + 1 - posS) / ml;
        PosM = posS - 1 + round(NR/2) * ml;
        
        [NPCR,Qmaxi] = FindPCRerr(NR,1,ml);
        
        for k = [-NPCR:-1,1:NPCR]
            w = 1 - abs(k) * Qmaxi;
            if w <= 0
                continue
            end
            AFASTAk = AFASTA{i}([1:PosM-1+k*ml,PosM:end]);
            ii = strcmp(AFASTA,AFASTAk);
            if any(ii)
                AMeshW(ii,i) = w;
            end
        end
    end
end
AMeshW = max(AMeshW,AMeshW.');