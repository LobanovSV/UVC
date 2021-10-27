function [FPOS,RS,RE,NR,isBad] = FindRegions(Allele,DistTol,AFTol, ...
    FASTA,posW,RAM,FILE)

% Real position R_real is pos0 + pos(1) - 2 + R
% Also, NV is (RS+1):(RE-1)

%% Check

N = length(Allele);
if N == 0
    FPOS = cell(0,1);
    RS = zeros(0,1);
    RE = zeros(0,1);
    NR = 0;
    isBad = false(0,1);
    return
end

%% Find FASTA positons and Regions with deviations from the FASTA

FPOS = cell(N,1);
RPOS = cell(N,2);
RLL = cell(N,2);
for i = 1:N
    [FPOS{i},RPOS(i,:),RLL(i,:),N2] = FindFRPOS(Allele{i},RAM,FILE,i);
end

%% FASTA position count

NA = 2 * N;
RC = zeros(1,N2);
RL = zeros(1,N2);
for i = 1:NA
    RC(RPOS{i}) = RC(RPOS{i}) + 1;
    RL(RPOS{i}) = max(RL(RPOS{i}),RLL{i}.');
end

%% Find STRs

[~,~,STR_S,STR_E] = FindSTR(FASTA);
STR_S = 2 - posW(1) + STR_S;
STR_E = 2 - posW(1) + STR_E;

%% Find regions

RS = find(RC > AFTol*NA).';
RL = 1 + DistTol * RL(RS);
NR = length(RS);
if NR == 0
    RE = RS;
    isBad = ~cellfun(@isempty,RPOS);
    return
end

% Find pair-wise distances
DD = true(1,N2);
for i = 1:length(STR_S)
    DD(max(1,STR_S(i)):min(N2,STR_E(i))) = false;
end
D = true(NR);
for i = 1:NR-1
    for j = i+1:NR
        d = sum(DD(RS(i):RS(j))) - 1;
        D(i,j) = d <= RL(i) && d <= RL(j);
    end
end
D = D & D.';

iS = ones(size(RS));
iE = ones(size(RS));
i = 1;
jp = 1;
while iE(i) < NR
    j = find(any(D(jp:iE(i),:),1),1,'last');
    jp = iE(i) + 1;
    if j > iE(i)
        iE(i) = j;
    else
        i = i + 1;
        iS(i) = jp;
        iE(i) = jp;
    end
end
RE = RS(iE(1:i));
RS = RS(iS(1:i));

% RS must be less than RE
ii = RS == RE;
if any(ii)
    RS(ii) = [];
    RE(ii) = [];
end
NR = length(RS);
if NR == 0
    isBad = ~cellfun(@isempty,RPOS);
    return
end

%% Find alleles with rare regions

isBad = false(N,2);
RSt = RS.';
REt = RE.';
for i = 1:NA
    if ~isempty(RPOS{i})
        isBad(i) = any(~any(RPOS{i} >= RSt & RPOS{i} <= REt,2));
    end
end