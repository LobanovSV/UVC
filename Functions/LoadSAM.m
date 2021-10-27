function [p1,p2,ID,NSEQ,SEQ,QUAL,POSS,POSE,PAIR,LEFT] = LoadSAM(ID, ...
    DirLoad,i,QUALminus,QUALdivide,RAM)

%% Parameters

Template = [ ...
    '%s' ... 1 QNAME String Query template NAME
    '%u16' ... 2 FLAG Int [0; 2^16-1] bitwise FLAG
    '%*s' ... 3 RNAME String Reference sequence NAME
    '%f' ... 4 POS Int [0; 2^31-1] 1-based leftmost mapping POSition
    '%u8' ... 5 MAPQ Int [0; 2^8-1] MAPping Quality
    '%s' ... 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
    '%*s' ... 7 RNEXT String Reference name of the mate/next read
    '%*s' ... 8 PNEXT Int [0; 2^31-1] Position of the mate/next read
    '%*s' ... 9 TLEN Int [-2^31 + 1; 2^31-1] observed Template LENgth
    '%s' ... 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
    '%s' ... 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
    '%*[^\n]'];

QUALminus = 33 + QUALminus;

NNTol = 0.1; % Tolerance for the number of N

%% Load

FileLoad = fullfile(DirLoad,[int2str(i) '.sam']);
File = fopen(FileLoad,'r');
s = textscan(File,Template,'Delimiter','\t');
fclose(File);

% Check second file
FileLoad = fullfile(DirLoad,[int2str(i) 'b.sam']);
if exist(FileLoad,'file')
    File = fopen(FileLoad,'r');
    sa = textscan(File,Template,'Delimiter','\t');
    fclose(File);
    for j = 1:length(s)
        s{j} = [s{j};sa{j}];
    end
    clear sa
end

%% Remove reads with zero mapping quality or reads with too many N

iiRem = s{4} == 0;
for j = 1:length(iiRem)
    iiRem(j) = iiRem(j) || mean(s{7}{j} <= QUALminus) > NNTol;
end
NRem = sum(iiRem); % Number of removed reads
if NRem ~= 0
    for j = 1:length(s)
        s{j}(iiRem) = [];
    end
end

%% Remove reads with more than 2 in a pair

[~,ii,PAIR] = unique(s{1},'stable');
sPAIR = sort(PAIR);
if length(sPAIR) > 2 && any(sPAIR(1:end-2)==sPAIR(3:end))
    jj = sPAIR(1:end-2)==sPAIR(3:end);
    sPAIR = unique(sPAIR(jj));
    iiRem = false(size(PAIR));
    for j = 1:length(sPAIR)
        iiRem = iiRem | PAIR == sPAIR(j);
    end
    for j = 1:length(s)
        s{j}(iiRem) = [];
    end
    NRem = NRem + sum(iiRem);
    [~,ii,PAIR] = unique(s{1},'stable');
end
LEFT = bitget(s{2},5) == 0;

%% Remove reads with the same edges

L1 = LEFT(ii);
L2 = ~L1;
jj = true(size(PAIR));
jj(ii) = false;
L2(PAIR(jj)) = LEFT(jj);
if any(L2==L1)
    sPAIR = find(L2==L1);
    iiRem = false(size(PAIR));
    for j = 1:length(sPAIR)
        iiRem = iiRem | PAIR == sPAIR(j);
    end
    for j = 1:length(s)
        s{j}(iiRem) = [];
    end
    NRem = NRem + sum(iiRem);
    [~,~,PAIR] = unique(s{1},'stable');
    LEFT = bitget(s{2},5) == 0;
end

%% Save number of removed reads

if NRem ~= 0
    ID = [ID ', ' int2str(NRem) ' reads were removed'];
end

%% Save data in the cell

NSEQ = length(s{1});
POSS = s{3};
SEQ = s{6};
QUAL = s{7};
POSE = zeros(NSEQ,1);
for j = 1:NSEQ
    POSE(j) = POSS(j) - 1 + length(SEQ{j});
    QUAL{j} = max(0,QUAL{j} - QUALminus) / QUALdivide;
end

%% Edges of all reads

if NSEQ == 0
    p1 = NaN;
    p2 = NaN;
else
    p1 = min(POSS);
    p2 = max(POSE);
end

%% Save

if ~RAM
    FileSave = fullfile(DirLoad,[int2str(i) '.mat']);
    save(FileSave,'SEQ','QUAL','PAIR','LEFT','POSS','POSE')
    SEQ = [];
    QUAL = [];
    POSS = [];
    POSE = [];
    PAIR = [];
    LEFT = [];
end