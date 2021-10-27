function [AlgD,QD,jjD,REF] = FindDiff(REF,FASTA,SEQ,QUAL,dEp,POSp, ...
    isPCRp,Ind,isN,IndN,posW,MinReads,Qmax,FracReads,TextType)

%% Some variables

NTT = length(TextType);
N2 = length(REF);
Np = length(dEp);

%% Find sequences in pairs

pos12 = zeros(2,2);
seq12 = cell(1,2);
qual12 = cell(1,2);
POSP = zeros(Np,2);
SEQP = cell(Np,1);
QUALP = cell(Np,1);
for i = 1:Np
    k = 1 + (dEp(i) < 0);
    % Skip PCRs dublications
    if isPCRp(i,k) || (~dEp(i) && any(isPCRp(i,:)))
        continue
    end
    for j = 1:2 % Left, right
        if ~Ind(i,j) % No left/right read
            pos12(j,:) = 0;
            continue
        end
        jj = 2*j-1:2*j; % Left or right
        
        % First (or the only one) allele
        pos = IndN{k}(POSp(i,jj,k));
        pp = pos(1):pos(2);
        isn = isN{k}(pp);
        if all(isn)
            seq = SEQ{Ind(i,j)};
            qual = QUAL{Ind(i,j)};
        else
            Nseq = 1 + pos(2) - pos(1);
            seq = repmat('*',1,Nseq);
            seq(isn) = SEQ{Ind(i,j)};
            qual = repmat(Qmax,1,Nseq);
            qual(isn) = QUAL{Ind(i,j)};
        end
        
        if ~dEp(i) % Two alleles
            pos2 = IndN{2}(POSp(i,jj,2));
            isn2 = isN{2}(pp);
            if pos(1)~= pos2(1) || pos(2) ~= pos2(2) || ~isequal(isn,isn2)
                isn2(pp<pos2(1)) = false;
                isn2(pp>pos2(2)) = false;
                i1 = 1 + sum(isN{2}(pos2(1):pos(1)-1));
                i2 = i1 - 1 + sum(isn2);
                seq2 = blanks(length(seq));
                seq2(pp>pos2(1) & pp<pos2(2)) = '*';
                seq2(isn2) = SEQ{Ind(i,j)}(i1:i2);
%                 disp([seq;seq2])
                col = find(seq==seq2);
                if isempty(col)
                    pos12(j,:) = 0;
                    break
                end
                ii = diff(col) ~= 1;
                colL = col([true,ii]);
                colR = col([ii,true]);
                [~,l] = max(colR + 1 - colL);
                i1 = colL(l);
                i2 = colR(l);
                
                pos = pos(1) - 1 + [i1,i2];
                seq = seq(i1:i2);
                qual = qual(i1:i2);
            end
        end
        
        pos12(j,:) = pos;
        seq12{j} = seq;
        qual12{j} = qual;
    end
    
    if all(pos12(:,1))
        pp = [min(pos12(:,1)),max(pos12(:,2))];
        Npp = pp(2) - pp(1) + 1;
        seq = blanks(Npp);
        qual = zeros(1,Npp);
        ii = pos12(1,1)+1-pp(1):pos12(1,2)+1-pp(1);
        seq(ii) = seq12{1};
        qual(ii) = qual12{1};
        
        seq2 = blanks(Npp);
        qual2 = zeros(1,Npp);
        ii = pos12(2,1)+1-pp(1):pos12(2,2)+1-pp(1);
        seq2(ii) = seq12{2};
        qual2(ii) = qual12{2};
        
        qual = max(qual,qual2);
        ii = seq ~= ' ' & seq2 ~= ' ' & seq ~= seq2;
        jj = seq == ' ';
        seq(jj) = seq2(jj);
        
        seq(ii) = ' ';
        qual(ii) = 0;
    else
        j = find(pos12(:,1),1);
        if isempty(j)
            continue
        end
        pp = [pos12(j,1),pos12(j,2)];
        seq = seq12{j};
        qual = qual12{j};
    end
    
    POSP(i,:) = pp;
    SEQP{i} = seq;
    QUALP{i} = qual;
end

%% Find number of nucleotides

Diff = zeros(NTT,N2);
% At least MinReads reads must have different starting/ending positions to
% be considered as non-dublicates
iiNum = zeros(NTT,N2);
ppS = zeros(NTT,N2,MinReads-1);
ppE = zeros(NTT,N2,MinReads-1);
for i = 1:Np
    if isempty(SEQP{i})
        continue
    end
    pp = POSP(i,1):POSP(i,2);
    
    ii = SEQP{i} == TextType;
    Diff(:,pp) = Diff(:,pp) + QUALP{i} .* ii;
    
    % Modify nucleotides with k-1 reads
    for k = 1:MinReads
        iik = find(ii & iiNum(:,pp)==(k-1)) + NTT * (pp(1)-1);
        if ~isempty(iik)
            for l = 1:k-1
                iil = iik + NTT * N2 * (l-1);
                iik = iik(ppS(iil) ~= pp(1) | ppE(iil) ~= pp(end));
                if isempty(iik)
                    break
                end
            end
            if ~isempty(iik)
                iiNum(iik) = iiNum(iik) + 1;
                if k ~= MinReads
                    iil = iik + NTT * N2 * (k-1);
                    ppS(iil) = pp(1);
                    ppE(iil) = pp(end);
                end
            end
        end
    end
end

%% Remove rare mutations

Diff(Diff<MinReads*Qmax) = 0; % Remove single mutations

Diff(Diff<FracReads*max(Diff,[],1)) = 0; % Remove single mutations

ii = 1:N2;
Diff(:,ii<posW(1)|ii>posW(2)) = 0; % Remove outer region

%% Replace nucleotides by asterics

Aster = FASTA == '*' & all(Diff(2:end,:)==0,1);

%% Find different nucleotides

TextType = TextType.';
jjD = false(1,N2);
ii = 0:NTT:NTT*(N2-1);
for i = 1:2
    [MaxDiff,iDiff] = max(Diff,[],1);
    Diff(iDiff+ii) = 0;
    REF1 = TextType(iDiff);
    jjD = jjD | (MaxDiff ~= 0 & REF1~=REF);
end
if any(Aster)
    jjD(Aster) = false;
    REF(Aster) = '*';
end
jjD = find(jjD);

%% Construct AlgD

N2D = length(jjD);
AlgD = repmat(blanks(N2D),Np,1);
QAlgD = zeros(Np,N2D);
for i = 1:Np
    if isempty(SEQP{i})
        continue
    end
    
    jj = jjD>=POSP(i,1) & jjD<=POSP(i,2);
    if any(jj)
        ii = 1-POSP(i,1)+jjD(jj);
        AlgD(i,jj) = SEQP{i}(ii);
        QAlgD(i,jj) = QUALP{i}(ii);
    end
end

% Remove alleles with unknown elements
% ii = any(AlgD(:)=='?'|AlgD(:)>='a',2);
% AlgD(ii) = ' ';
% QAlgD(ii) = 0;

% Find unique alleles
[AlgD,~,iiD] = unique(AlgD,'rows');
N1 = size(AlgD,1);
QD = zeros(N1,N2D); % Occurence
for i = 1:N1
    ii = iiD==i;
    QD(i,:) = sum(QAlgD(ii,:),1);
end

% Remove empty line
if all(AlgD(1,:) == ' ')
    AlgD(1,:) = [];
    QD(1,:) = [];
    if isempty(AlgD)
        jjD = [];
    end
end

% wS = zeros(size(ND));
% wN = zeros(size(ND));
% for i = 1:size(AlgD,2)
%     for j = 2:length(TextType)
%         if Diff(j,i) ~= 0
%             ii = AlgD(:,i) == TextType(j);
%             wN(ii) = wN(ii) + 1;
%             wS(ii) = wS(ii) + Diff(j,i);
%         end
%     end
% end
% ND = ND .* wS ./ wN;
% ND(isnan(ND)) = 0;
end