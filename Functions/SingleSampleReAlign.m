function [Allele,Alg,Err] = SingleSampleReAlign(g,rf,SEQ,QUAL, ...
    POSS,POSE,PAIR,LEFT,posW,dPos,dPOS,RAM,FILE,i,AlleleOnly)

%% Parameters

FracReads = eps;
TextType = '*ACGT'.';

%% Load data

if ~RAM
    File = fullfile(DirPCorUNIX,'BAMs','Temp',FILE,[int2str(i) '.mat']);
    load(File,'SEQ','QUAL','PAIR','LEFT','POSS','POSE')
end

%% Best quality

Qmax = max(cellfun(@max,QUAL));
maxSL = 1 + max(POSE - POSS);

%% Find Alleles and REF

% A1
Allele = cell(4,1);
Allele{1} = strjoin(rf(1,:),'');

% FASTA
FASTA = strjoin(rf(end,:),'');
posFASTA = cumsum(FASTA~='*');
posW(1) = find(posFASTA==posW(1),1);
posW(2) = find(posFASTA==posW(2),1);

% REF
REF = FASTA;
ii = REF == '*';
REF(ii) = Allele{1}(ii);

if all(g(1,1,:)==g(1,2,:))
    Allele{2} = '';
else
    Allele{2} = strjoin(rf(2,:),'');
    ii = REF == '*';
    REF(ii) = Allele{2}(ii);
end
Allele{3} = REF;
Allele{4} = FASTA;

%% Find indexes of pairs

Np = max(PAIR);
Ind = zeros(Np,2);
Indt = 1:length(PAIR);
Ind(PAIR(LEFT),1) = Indt(LEFT);
Ind(PAIR(~LEFT),2) = Indt(~LEFT);

%% Create Final Alignment and find error

[Err,Alg,Allele] = CreateFinalAlignment(Allele,SEQ,QUAL,POSS,PAIR,LEFT, ...
    Ind,posW,maxSL,dPos,dPOS,FracReads,Qmax,TextType);

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