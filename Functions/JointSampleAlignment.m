function [RS,RE,A,AisN,G,Gp,IDL,IDLp,Allele,Alg,Err] = ...
    JointSampleAlignment(Allele,Alg,Err,SEQ,QUAL,POSS,POSE,PAIR, ...
    LEFT,FASTA,posW,dPos,dPOS,AFTol,Fast,RAM,FILE,AlleleOnly)

% IDL contains InDel lengths shifted to the reference value with the
% following order: Sum, Max, Min, Diff, Add3, Dom3, Rec3, Het3

%% Default parameters

% Tollerance for region frequency
if nargin < 14 || isempty(AFTol)
    AFTol = 1e-3;
end
if nargin < 15 || isempty(Fast)
    Fast = false;
end
if nargin < 16 || isempty(RAM)
    RAM = true;
    FILE = '';
    AlleleOnly = false;
end

%% Parameters

DistTol = 2; % Tollerance for the distance between regions
ALTol = 0.5; % Must be at least one allele with GL above the threshold

%% Find regions of interest

[FPOS,RS,RE,NR,isBad] = FindRegions(Allele,DistTol,AFTol,FASTA,posW, ...
    RAM,FILE);

%% Find alleles and genotypes

MergeR = true;
while MergeR
    MergeR = false;
    N = length(Allele);
    G = zeros(N,2,NR);
    A = cell(NR,1);
    AisN = cell(NR,1);
    AMeshW = cell(NR,1);
    iRemSingle = false(NR,1);
    RF = FindRegions_RF(RS,RE,NR,FASTA,posW);
    for i = 1:NR
        [A{i},AisN{i},G(:,:,i),AMeshW{i}] = FindAandG(RS(i),RE(i), ...
            Allele,FPOS,isBad,FASTA,posW,RAM,FILE);
        iRemSingle(i) = length(A{i}) == 1;
        
        % Check whether regions should be merged
        if i ~= 1
            n1 = length(A{i-1});
            n2 = length(A{i});
            S = cell(n1,n2);
            for i1 = 1:n1
                for i2 = 1:n2
                    S{i1,i2} = [A{i-1}{i1},RF{2*i-1},A{i}{i2}];
                end
            end
            if length(unique(S(:))) ~= n1*n2
                RS(i) = [];
                RE(i-1) = [];
                NR = NR - 1;
                MergeR = true;
                break
            end
        end
    end
end
isBad = any(isBad,2);

% Remove single alleles
if any(iRemSingle)
    RS(iRemSingle) = [];
    RE(iRemSingle) = [];
    G(:,:,iRemSingle) = [];
    A(iRemSingle) = [];
    AisN(iRemSingle) = [];
    AMeshW(iRemSingle) = [];
    NR = length(RS);
end
Gp = zeros(N,NR);

if NR == 0
    IDL = [];
    IDLp = [];
    return
end

clear n1 n2 S i1 i2 MergeR

%% Find probabilities and correct genotypes

for iCorrG = 1:10
    Break = true;
    RF = FindRegions_RF(RS,RE,NR,FASTA,posW);
    rf = repmat(RF,3,1); % A1, A2, REF, FASTA
    
    IC = cell(NR,1);
    for j = 1:NR
        IC{j} = eye(length(A{j}),1);
    end
    
    isGood = ~isBad;
    for i = 1:N
        clc
        fprintf(['Find probabilities and correct genotypes: ' ...
            'iteration %u, sample %u from %u\n'],iCorrG,i,N)
        
        g = G(i,:,:);
        [g,Gp(i,:),Ap] = FindGp(g,A,RS,RE,RF,SEQ{i},QUAL{i},POSS{i}, ...
            POSE{i},PAIR{i},LEFT{i},dPos,dPOS,RAM,FILE,i);
        for j = 1:NR
            IC{j} = IC{j} | AMeshW{j} * Ap{j} > ALTol;
        end
        if Fast
            G(i,:,:) = g;
            continue
        end
        if isBad(i) || ~isequal(g,G(i,:,:))
            for j = 1:NR
                j2 = 2*j;
                g1 = g(1,1,j);
                g2 = g(1,2,j);
                if (g1 == 1 && g2 == 1) || (length(A{j}{1}) == 1 && ...
                        length(A{j}{g1}) == 1 && length(A{j}{g2}) == 1)
                    rf{1,j2} = A{j}{g1};
                    rf{2,j2} = A{j}{g2};
                    rf{3,j2} = A{j}{1};
                else
                    if isGood(i) && isequal(G(i,:,j),g(1,:,j))
                        k = i;
                    else
                        k = find(isGood & all(G(:,:,j)==g(1,:,j),2),1);
                    end
                    if isempty(k)
                        IDL_Ind = [g1,g2,1];
                        aisn = AisN{j}(IDL_Ind,:);
                        aisn0 = cumsum(max(aisn));
                        for m = 1:3
                            isN = false(1,aisn0(end));
                            for l = 1:length(aisn0)
                                lm2 = aisn0(l);
                                lm1 = lm2 + 1 - aisn(m,l);
                                if lm1 <= lm2
                                    isN(lm1:lm2) = true;
                                end
                            end
                            rf{m,j2} = repmat('*',size(isN));
                            rf{m,j2}(isN) = A{j}{IDL_Ind(m)};
                        end
                    else
                        [Allelek,FPOSk] = LoadAlleleFPOS(Allele, ...
                            FPOS,RAM,FILE,k);
                        ii = FPOSk(RS(j))+1:FPOSk(RE(j))-1;
                        rf{1,j2} = Allelek{1}(ii);
                        if isempty(Allelek{2})
                            rf{2,j2} = Allelek{1}(ii);
                        else
                            rf{2,j2} = Allelek{2}(ii);
                        end
                        rf{3,j2} = Allelek{end}(ii);
                    end
                end
            end
            % Re-align
            [Allele{i},Alg{i},Err(i)] = SingleSampleReAlign(g, ...
                rf,SEQ{i},QUAL{i},POSS{i},POSE{i},PAIR{i},LEFT{i}, ...
                posW,dPos,dPOS,RAM,FILE,i,AlleleOnly);
            % Find FASTA positons
            FPOS{i} = FindFRPOS(Allele{i},RAM,FILE,i);
            % Save genotype
            G(i,:,:) = g;
        end
    end
    isBad(:) = false;
    
    %% Remove zero-called alleles
    
    for i = NR:-1:1
        NAi = length(A{i});
        
        % Replace non-reliable alleles by nearest
        for j = find(~IC{i}).'
            ii = any(G(:,:,i)==j,2);
            if any(ii)
                G(G(:,1,i)==j,1,i) = 1;
                G(G(:,2,i)==j,2,i) = 1;
                isBad(ii) = true;
            end
        end
        
        % Remove zero-called alleles
        g = reshape(G(:,:,i),2*N,1);
        [gu,~,g] = unique([g;1]);
        if length(gu) == 1
            RS(i) = [];
            RE(i) = [];
            A(i) = [];
            AisN(i) = [];
            AMeshW(i) = [];
            G(:,:,i) = [];
            Gp(:,i) = [];
            Break = false;
        elseif length(gu) < NAi
            A{i} = A{i}(gu);
            AisN{i} = AisN{i}(gu,:);
            AMeshW{i} = AMeshW{i}(gu,gu);
            G(:,:,i) = reshape(g(1:end-1),N,2);
            [A{i},AisN{i},RS(i),RE(i)] = RemoveTheSameSeq(A{i},AisN{i}, ...
                RS(i),RE(i));
            Break = false;
        end
    end
    NR = length(RS);
    if NR == 0
        G = zeros(N,2,0);
        Gp = zeros(N,0);
        IDL = [];
        IDLp = [];
        return
    end
    if Fast || Break
        RF = FindRegions_RF(RS,RE,NR,FASTA,posW);
        break
    end
end

%% Find IDL (InDel lengths shifted to the reference value)

% Find lengths
[IDL_L,IDL_Ind] = Find_IDL_Ind(A);

IDL = zeros(N,8,NR);
IDLp = zeros(N,8,NR);
for i = 1:N
    clc
    fprintf(['Find InDel lengths shifted to the reference value: ' ...
        'sample %u from %u\n'],i,N)
    [IDL(i,:,:),IDLp(i,:,:)] = FindIDLp(G(i,:,:),A,IDL_L,IDL_Ind,RF, ...
        SEQ{i},QUAL{i},POSS{i},POSE{i},PAIR{i},LEFT{i},dPos,dPOS,RAM, ...
        FILE,i);
end