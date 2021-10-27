function [pos0,chr,posW,FASTA,ID,NSEQ,SEQ,QUAL,POSS,POSE,PAIR,LEFT] = ...
    LoadSAMs(DirLoad,ID,dPOS,QUALminus,QUALdivide,RAM,RefGenome)

%% Default variables

if nargin < 2 || isempty(dPOS)
    dPOS = 200;
end
if nargin < 4 || isempty(QUALminus)
    QUALminus = 10;
end
if nargin < 5 || isempty(QUALdivide)
    QUALdivide = 2;
end
if nargin < 6 || isempty(RAM)
    RAM = true;
end
if nargin < 7 || isempty(RefGenome)
    RefGenome = 'hg19';
end

%% DirLoad

FileLoad = fullfile(DirLoad,'Pos.mat');
load(FileLoad,'chr','posW')

%% Load data

N = length(ID);
p1 = zeros(N,1);
p2 = zeros(N,1);
NSEQ = zeros(N,1);
SEQ = cell(N,1);
QUAL = cell(N,1);
POSS = cell(N,1);
POSE = cell(N,1);
PAIR = cell(N,1);
LEFT = cell(N,1);
for i = 1:N
    [p1(i),p2(i),ID{i},NSEQ(i),SEQ{i},QUAL{i},POSS{i},POSE{i}, ...
        PAIR{i},LEFT{i}] = LoadSAM(ID{i},DirLoad,i,QUALminus, ...
        QUALdivide,RAM);
end
p1 = nanmin(p1);
p2 = nanmax(p2);

%% Shift and truncate further

% Load FASTA
if ~ischar(chr)
    chr = int2str(chr);
end
load(fullfile(cd,'RefGenome',RefGenome,[chr '.mat']),'FASTA')

p1 = min(posW(1),p1 - dPOS);
p2 = max(posW(2),p2 + dPOS);
FASTA = FASTA(p1:p2);
pos0 = p1 - 1;
posW = posW - pos0;
for i = 1:N
    [POSS{i},POSE{i}] = SubtractPOS(POSS{i},POSE{i},pos0,RAM,DirLoad,i);
end
end

function [POSS,POSE] = SubtractPOS(POSS,POSE,pos0,RAM,DirLoad,i)
% Load data
if ~RAM
    File = fullfile(DirLoad,[int2str(i) '.mat']);
    load(File,'POSS','POSE')
end

% Subtract
POSS = POSS - pos0;
POSE = POSE - pos0;

% Save data
if ~RAM
    save(File,'POSS','POSE','-append')
    POSS = [];
    POSE = [];
end
end