function [row,colL,colR,S] = Align2Right(L0Del,L0Ins,row0,colL0,colR0, ...
    M,Row,ColL,ColR,LL,reverse,SThr)
%% Additional variables

Ns = size(M,2);

%% Default values and checks

row = [];
colL = [];
colR = [];
S = 0;
% End is finished
if (reverse && colL0 == 1) || (~reverse && colR0 == Ns)
    return
end

%% Truncate

NL = length(LL);
for i = 1:NL
    if reverse
        ii = ColL{i} <= colL0;
    else
        ii = ColR{i} >= colR0;
    end
    if (reverse && isinf(L0Ins)) || (~reverse && isinf(L0Del))
        ii = ii & Row{i} > row0;
    end
    if (~reverse && isinf(L0Ins)) || (reverse && isinf(L0Del))
        ii = ii & Row{i} < row0;
    end
    Row{i} = Row{i}(ii);
    ColL{i} = ColL{i}(ii);
    ColR{i} = ColR{i}(ii);
    LL{i} = LL{i}(ii);
end

%% Reverse

if reverse
    M = M(:,end:-1:1);
    
    colR0 = Ns + 1 - colL0;
    
    for i = 1:NL
        Temp = ColL{i};
        ColL{i} = Ns + 1 - ColR{i};
        ColR{i} = Ns + 1 - Temp;
    end
    
    Temp = L0Ins;
    L0Ins = L0Del;
    L0Del = Temp;
end

%% Align

m = 0;
MaxStep = zeros(1,1);
MStep = zeros(1,0); % Mismatch step: odd - 0m,1m,2m; even - +0m,+1m,...
Smax = min(SThr*Ns,sum(~M(row0,colR0+1:end))); % Maximum mismatch
POS = zeros(0,2);
SCO = zeros(0,1);
while m < Smax
    EpsTol = 10 * (m + 1) * eps;
    for i = 1:length(MaxStep)
        [row,colL,colR,MaxStep(i),POS,SCO,S] = SinglePiece(L0Del,L0Ins, ...
            -m-1,row0,colR0,Row,ColL,ColR,LL,Ns,POS,SCO,0,MStep(i,:),M, ...
            reverse,EpsTol);
        % Break the loop if alignment is found
        if ~isempty(row)
            break
        end
    end
    if ~isempty(row)
        if length(S) ~= 1 || S < -1
            error('Wrong S')
        end
        S = S + m + 1;
        break
    end
    if m == 0
        MStep = (2:2*MaxStep).';
    else
        MStep = repelem(MStep,max(1,2*MaxStep+1-MStep(:,end)),1);
        MStep = [MStep zeros(size(MStep,1),1)]; %#ok
        mstep = [];
        iRem = false(size(MStep,1),1);
        for i = 1:size(MStep,1)
            if ~isequal(mstep,MStep(i,:))
                mstep = MStep(i,:);
                j = MStep(i,end-1);
                if mod(j,2) == 1 && m>NL-2 && ...
                        all(MStep(i,end+1-NL:end-1)==j)
                    iRem(i) = true;
                end
            else
                j = j + 1;
            end
            MStep(i,end) = j;
        end
        if any(iRem)
            MStep(iRem,:) = [];
        end
    end
    MaxStep = zeros(size(MStep,1),1);
    m = m + 1;
end

%% No alignment was found

if isempty(row)
    S = sum(~M(row0,colR0+1:end));
    row = row0;
    colL = colR0+1;
    colR = Ns;
elseif colL(1) > colR0+1 % The first element might be shifted
    row = [row0; row];
    colR = [colL(1)-1; colR];
    colL = [colR0+1; colL];
end

%% Reverse before return

if reverse
    Temp = colL;
    colL = Ns + 1 - colR;
    colR = Ns + 1 - Temp;
    
    row = row(end:-1:1);
    colL = colL(end:-1:1);
    colR = colR(end:-1:1);
end