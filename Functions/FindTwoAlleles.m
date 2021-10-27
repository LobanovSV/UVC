function AlleleD = FindTwoAlleles(AlgD,QD,TextType,REF,Qmax)

NU = size(AlgD,1);
AlleleD = {REF;REF};
if NU == 1 && all(AlgD==' ')
    return
end

%% Find weight Q

Q = zeros(NU);
for i = 1:NU
    for j = 1:i-1
        NE = AlgD(i,:) ~= AlgD(j,:);
        Q(i,j) = sum(QD(i,NE) .* QD(j,NE));
    end
end
Q = Q + Q.';

%% Find indeces

i0 = 1:NU;
i1 = [];
i2 = [];
while ~isempty(i0)
    w = sum(Q(i1,i0),1)-sum(Q(i2,i0),1);
    [~,j] = max(abs(w));
    if w(j) == 0
        Q00 = Q(i0,i0);
        [Q00Max,j] = max(Q00(:));
        if Q00Max == 0
            i1 = [i1;i0.']; %#ok
            i2 = [i2;i0.']; %#ok
            break
        end
        [i,j] = ind2sub(size(Q00),j);
        i1 = [i1;i0(i)]; %#ok
        i2 = [i2;i0(j)]; %#ok
        i0([i,j]) = [];
        continue
    elseif w(j) > 0
        i2 = [i2;i0(j)]; %#ok
    else
        i1 = [i1;i0(j)]; %#ok
    end
    i0(j) = [];
end

%% Find two alleles

NTT = length(TextType);
N2 = size(AlgD,2);
Diff1 = zeros(NTT,N2);
Diff2 = zeros(NTT,N2);
for i = 1:NTT
    Diff1(i,:) = sum((AlgD(i1,:)==TextType(i)) .* QD(i1,:),1);
    Diff2(i,:) = sum((AlgD(i2,:)==TextType(i)) .* QD(i2,:),1);
end

[MaxDiff,iDiff1] = max(Diff1,[],1);
iDiff1(MaxDiff<Qmax) = 0;
[MaxDiff,iDiff2] = max(Diff2,[],1);
iDiff2(MaxDiff<Qmax) = 0;

for i = 1:NTT
    ii = iDiff1 == i;
    AlleleD{1}(ii) = TextType(i);
    ii = iDiff2 == i;
    AlleleD{2}(ii) = TextType(i);
end

if sum(AlleleD{1}=='*') > sum(AlleleD{2}=='*')
    AlleleD = AlleleD([2 1]);
end
end