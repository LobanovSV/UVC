function [PCR,isN,IndN] = FindPCRsREF(Allele,STR,Qmax)

isN = cell(2,1);
IndN = cell(2,1);
PCR = struct;
for i = 1:2
    if isempty(Allele{i})
        continue
    end
    
    NSTR = length(STR(i).PosS);
    N = 1 + 2 * sum(STR(i).NPCR);
    PCR(i).seq = cell(N,1);
    PCR(i).pos = STR(i).PosS - 1;
    PCR(i).dE = zeros(N,1);
    PCR(i).dN = zeros(NSTR,N);
    PCR(i).Ind = zeros(N,1);
    IndN{i} = cell(N,1);
    
    j = 1;
    isN{i} = repmat(Allele{i} ~= '*',N,1);
    PCR(i).seq{1} = Allele{i}(isN{i}(1,:));
    IndN{i}{1} = find(isN{i}(1,:));
    
    for k = 1:NSTR
        ml = STR(i).ML(k);
        poss = STR(i).PosS(k);
        pose = STR(i).PosE(k);
        possReal = IndN{i}{1}(poss);
        [npcr,Qmaxi] = FindPCRerr((pose + 1 - poss) ./ ml,Qmax,ml);
        for l = [-npcr:-1,1:npcr]
            j = j + 1;
            PCR(i).dE(j) = abs(l) * Qmaxi;
            PCR(i).seq{j} = PCR(i).seq{1}([1:poss-1+ml*(npcr-l), ...
                poss+ml*npcr:end]);
            if l < 0 % Insertion
                isN{i}(j,possReal+l*ml:possReal-1) = true;
            else % Deletion
                isN{i}(j,possReal:IndN{i}{1}(poss-1+l*ml)) = false;
            end
            IndN{i}{j} = find(isN{i}(j,:));
            PCR(i).dN(k,j) = -l*ml;
            PCR(i).Ind(j) = k;
        end
    end
end