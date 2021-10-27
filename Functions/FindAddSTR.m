function [Allele,posW,STR] = FindAddSTR(Allele,posW)

%% Parameters

%% Find STR

STR = struct;
for i = 1:2
    if isempty(Allele{i})
        continue
    end
    
    isn = Allele{i} ~= '*';
    indn = find(isn);
    
    % Find STRs
    post = [find(indn<=posW(1),1,'last'),find(indn>=posW(2),1,'first')];
    [Motif,ML,PosS,PosE] = FindSTR(Allele{i}(isn),post);
    
    % Number of PCRs motifs
    NPCR = FindPCRerr((PosE + 1 - PosS) ./ ML);
    
    % Write to STR
    STR(i).Motif = Motif;
    STR(i).ML = ML;
    STR(i).PosS = PosS;
    STR(i).PosE = PosE;
    STR(i).NPCR = NPCR;
    if isempty(ML)
        continue
    end
    
    % Expand
    for j = 1:length(ML)
        poss = indn(PosS(j));
        if PosS(j) == 1
            posp = 0;
        else
            posp = indn(PosS(j)-1);
        end
        nadd = NPCR(j) * ML(j) - (poss - posp - 1);
        if nadd > 0
            for k = 1:length(Allele)
                if ~isempty(Allele{k})
                    Allele{k} = [Allele{k}(1:posp), ...
                        repmat('*',1,nadd),Allele{k}(posp+1:end)];
                end
            end
            ii = posW > posp;
            posW(ii) = posW(ii) + nadd;
            isn = [isn(1:posp),false(1,nadd),isn(posp+1:end)];
            indn(PosS(j):end) = indn(PosS(j):end) + nadd;
        end
    end
end