function [row,colL,colR,sco] = AlignReadLMR(L0Del,L0Ins,M,Row,ColL, ...
    ColR,LL,k)
% Align Left-Middle-Right part of the read

%% Parameters

SThr = 0.1;

%% Find the longest piece i0,[j0l,j0r]

row0 = Row{1}(k);
colL0 = ColL{1}(k);
colR0 = ColR{1}(k);

%% Align to right

% Starting position
[rowr,colLr,colRr,scor] = Align2Right(L0Del,L0Ins,row0,colL0,colR0,M, ...
    Row,ColL,ColR,LL,false,SThr);

%% Align to left

% Starting position
[rowl,colLl,colRl,scol] = Align2Right(L0Del,L0Ins,row0,colL0,colR0,M, ...
    Row,ColL,ColR,LL,true,SThr);

%% Matching quality

sco = scol + scor;

%% Combine

% Remove pieces with poor quality - right
for i = length(rowr):-1:1
    if mean(~M(rowr(i),colLr(i):colRr(i))) >= SThr
        if i == 1
            rowr = row0;
        else
            colRr(i-1) = colRr(i);
            rowr(i) = [];
            colLr(i) = [];
            colRr(i) = [];
        end
    else
        break
    end
end
% Merge with right
if ~isempty(rowr) && rowr(1) == row0
    colR0 = colRr(1);
    rowr(1) = [];
    colLr(1) = [];
    colRr(1) = [];
end

% Remove pieces with poor quality - left
for i = 1:length(rowl)
    if mean(~M(rowl(1),colLl(1):colRl(1))) >= SThr
        if length(rowl) == 1
            rowl = row0;
        else
            colLl(2) = colLl(1);
            rowl(1) = [];
            colLl(1) = [];
            colRl(1) = [];
        end
    else
        break
    end
end
% Merge with right
if ~isempty(rowl) && rowl(end) == row0
    colL0 = colLl(end);
    rowl(end) = [];
    colLl(end) = [];
    colRl(end) = [];
end

% Merge
row = [rowl; row0; rowr];
colL = [colLl; colL0; colLr];
colR = [colRl; colR0; colRr];

% clearvars -except sco row colL colR M

%% Shift matches to left

i = length(row);
while i > 1
    if row(i) == row(i-1)
        colL(i) = colL(i-1);
        i = i - 1;
        row(i) = [];
        colL(i) = [];
        colR(i) = [];
    elseif colR(i-1) < colL(i-1)
        i = i - 1;
        row(i) = [];
        colL(i) = [];
        colR(i) = [];
    elseif M(row(i),colL(i)-1)
        colL(i) = colL(i) - 1;
        colR(i-1) = colR(i-1) - 1;
    else
        i = i - 1;
    end
end
if length(row) > 1 && row(1) == row(2)
    row(1) = [];
    colL(1) = [];
    colR(1) = [];
    colL(1) = 1;
end

%% Some checks

Ns = size(M,2);
if any(colR<colL)
    error('Consider')
end
if colL(1) ~= 1 || colR(end) ~= Ns || any(colR(1:end-1)+1 ~= colL(2:end))
    error('Consider')
end
if any(diff(row)==0)
    error('Consider')
end