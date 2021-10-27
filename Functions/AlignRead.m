function [posS,posE,sco,InDel] = AlignRead(seq,ML,pos0,L0Ins,L0Del)
% 1) pos = pos0 - i + j for M_{i,j}
%    ML - matching matrix with letters
%    M - matching matrix, i.e. ML==seq
% 2) Reads with low quality are not removed but must be aligned without
%    insertions or deletions
% Parameters

% Score S = \Delta i \exp^{L0-L}, L0 = 4 gives 1% error

%% Optional input variables

% Threshold for insertions should be greater than threshold for deletions
% because insertions modify REF
if nargin < 4 || isempty(L0Ins)
    L0Ins = 5;
end
% Threshold for deletions
if nargin < 5 || isempty(L0Del)
    L0Del = 5;
end

%% Parameters

% Discrepancy in LL0 to be checked
dLL0 = 3;

%% Calculate matching matrix

if any(seq=='*')
    seq(seq=='*') = [];
end
Ns = length(seq);
if Ns < size(ML,2)
    ML = ML(:,1:Ns);
elseif Ns > size(ML,2)
    error('Wrong maxSL!')
end
% SaveAlgCoins(seq,ML,'Read')
M = ML == seq;

%% If no insertion or deletions are permitted

if isinf(L0Ins) && isinf(L0Del)
    [sco,row] = min(sum(~M,2));
    posS = pos0 - row + 1;
    posE = pos0 - row + Ns;
    return
end

%% Find pieces

Row = cell(3,1);
ColL = cell(3,1);
ColR = cell(3,1);
LL = cell(3,1);

% 0 mismatch: L
[col,row] = find(M.');
ii = diff(row) ~= 0 | diff(col) ~= 1;
ColL{1} = col([true; ii]);
Row{1} = row([true; ii]);
ColR{1} = col([ii; true]);
LL{1} = ColR{1} + 1 - ColL{1};
clear row col

% 1 mismatch: L-m-L
ii = find(Row{1}(2:end)==Row{1}(1:end-1) & ...
    ColL{1}(2:end)-2==ColR{1}(1:end-1));
Row{2} = Row{1}(ii);
ColL{2} = ColL{1}(ii);
ColR{2} = ColR{1}(ii+1);
LL{2} = ColR{2} - ColL{2};

% 2 mismatch: L-mm-L
ii = find(Row{1}(2:end)==Row{1}(1:end-1) & ...
    ColL{1}(2:end)-3==ColR{1}(1:end-1));
Row{3} = Row{1}(ii);
ColL{3} = ColL{1}(ii);
ColR{3} = ColR{1}(ii+1);

% 2 mismatch: L-m-L-m-L
ii = find(Row{2}(2:end)==Row{2}(1:end-1) & ...
    ColL{2}(2:end)<=ColR{2}(1:end-1));
Row{3} = [Row{3}; Row{2}(ii)];
ColL{3} = [ColL{3}; ColL{2}(ii)];
ColR{3} = [ColR{3}; ColR{2}(ii+1)];
LL{3} = ColR{3} - ColL{3} - 1;

%% Truncate because pieces with length L<=L0 cannot be considered

for i = 1:length(LL)
    ii = LL{i} > min(L0Del,L0Ins);
    Row{i} = Row{i}(ii);
    ColL{i} = ColL{i}(ii);
    ColR{i} = ColR{i}(ii);
    LL{i} = LL{i}(ii);
end

%% Perform alignment

LLmax = max(LL{1});
sco = inf;
for k = find(LL{1}>=LLmax-dLL0).'
    [rowt,colLt,colRt,scot] = AlignReadLMR0(L0Del,L0Ins,M,Row,ColL, ...
        ColR,LL,k);
    if scot < sco
        row = rowt;
        colL = colLt;
        colR = colRt;
        sco = scot;
    end
end

%% Save alignment

posS = pos0 - row(1) + colL(1);
posE = pos0 - row(end) + colR(end);
if length(row) == 1
    InDel = zeros(0,2);
else
    InDel = pos0 + [colR(1:end-1)-row(1:end-1),colL(2:end)-row(2:end)];
end

pS = InDel(:,2);
pE = InDel(:,1);
ii = pE < pS;
pS(ii) = InDel(ii,1) + 1;
pE(ii) = InDel(ii,2) - 1;
for i = 1:length(pS)
    if sum(pE<=pS(i) & pS>=pE(i)) > 1
        FileSave = ['Error_' int2str(randi(1e7)) '.mat'];
        save(FileSave)
        error(['Data is saved in ' FileSave])
    end
end