function [Allele,Alg,Err] = FindAlignment(INDEL,InDel,REF,SEQ,QUAL, ...
    POSS,PAIR,LEFT,dEp,POSp,isPCRp,Ind,posW,maxSL,dPos,dPOS,Qmax, ...
    MinReads,FracReads)

%% No reads

if isempty(SEQ)
    Allele = {REF;'';'';''};
    Alg = repmat({repmat(blanks(length(REF)),0,1)},4,1);
    Err = struct('nE_mean',0,'NE_max',0,'nE_max',0,'sEwrong',0,'sE',0, ...
        'nSNP',0,'Qmax',0);
    return
end

%% Parameters

TextType = '*ACGT'.';

%% Convert INDEL to Del

[REF,FASTA,posW,isN,IndN] = INDEL2DEL(REF,posW,INDEL,InDel);

%% Do not consider homozygous InDels twice

HomoInDel = ~any(xor(InDel(:,1),InDel(:,2)));
if HomoInDel
    dEp(:) = inf;
end

%% Find different alleles

[AlgD,QD,jjD,REF] = FindDiff(REF,FASTA,SEQ,QUAL,dEp,POSp,isPCRp,Ind, ...
    isN,IndN,posW,MinReads,Qmax,FracReads,TextType);

Allele = {REF;'';REF;FASTA};
if ~isempty(jjD)
    % Find two alleles
    AlleleD = FindTwoAlleles(AlgD,QD,TextType,REF(jjD),MinReads*Qmax);
    
    Allele{1}(jjD) = AlleleD{1};
    if ~strcmp(AlleleD{1},AlleleD{2})
        Allele{2} = REF;
        Allele{2}(jjD) = AlleleD{2};
    end
end

%% Create Final Alignment and find error

[Err,Alg,Allele] = CreateFinalAlignment(Allele,SEQ,QUAL,POSS,PAIR,LEFT, ...
    Ind,posW,maxSL,dPos,dPOS,FracReads,Qmax,TextType);