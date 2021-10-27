function [REF,FASTA,posW,isN,IndN] = INDEL2DEL(REF,posW,INDEL,InDel)

%% Include InDels

NREF = length(REF);
Ref = cell(3,1);
Ind = cell(3,1);
I = zeros(3,1);
NStar = cell(3,1);
for j = 1:2
    [Ref{j},Ind{j}] = INDEL2DEL_1(REF,INDEL(InDel(:,j),:));
    I(j) = length(Ref{j});
    NStar{j} = zeros(1,I(j));
end
Ref{3} = REF;
Ind{3} = 1:NREF;
I(3) = NREF;
NStar{3} = zeros(1,NREF);

%% Find insertions/deletions

while all(I ~= 1)
    MaxInd = max([Ind{1}(I(1)),Ind{2}(I(2)),Ind{3}(I(3))]);
    for j = 1:3
        if Ind{j}(I(j)) == MaxInd
            I(j) = I(j) - 1;
        else
            NStar{j}(I(j)) = NStar{j}(I(j)) + 1;
        end
    end
end

%% Add deletions

for j = 1:3
    ii = find(NStar{j});
    for i = length(ii):-1:1
        Ref{j} = [Ref{j}(1:ii(i)),repmat('*',1,NStar{j}(ii(i))), ...
            Ref{j}(1+ii(i):end)];
    end
end

%% Find FASTA, REF, posW, isN

FASTA = Ref{3};
REF = FASTA;
iiFASTA = cumsum(FASTA ~= '*');
isN = cell(1,2);
IndN = cell(1,2);
for i = 1:2
    ii = REF == '*';
    REF(ii) = Ref{i}(ii);
    
    posW(i) = find(iiFASTA == posW(i),1);
    
    isN{i} = Ref{i} ~= '*';
    IndN{i} = find(isN{i});
end

%% Checks

if any(REF=='*') || length(Ref{1}) ~= length(Ref{2}) || ...
        length(Ref{1}) ~= length(REF)
    error('Wrong')
end

end

function [REF,Ind] = INDEL2DEL_1(REF,INDEL)

Ind = 1:length(REF);

if isempty(INDEL)
    return
end

%% Replace

N = size(INDEL,1);
for i = N:-1:1
    p1 = INDEL(i,1);
    p2 = INDEL(i,2);
    if p1 < p2
        REF = REF([1:p1,p2:end]);
        Ind = Ind([1:p1,p2:end]);
    else
        REF = REF([1:p1,p2:p1,p1+1:end]);
        Ind = Ind([1:p1,p2:p1,p1+1:end]);
    end
end
end