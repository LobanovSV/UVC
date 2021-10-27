function [row,colL,colR,Step,POS,SCO,S] = SinglePiece(L0Del,L0Ins,sco, ...
    row0,col0,Row,ColL,ColR,LL,Ns,POS,SCO,Step,MStep,M,reverse,EpsTol)
%% Parameters

%% Calculate step and current position

Step = Step + 1;
pos0 = [row0, col0];

%% Default values if the end is not found

row = [];
colL = [];
colR = [];
S = 0;

%% Shift starting position

if Step == 1
    N2 = sum(MStep == 2);
    if N2 ~= 0
        sco = sco + N2;
        kmis = find(~M(row0,col0+2:end),N2,'first');
        if length(kmis) ~= N2
            error('Impossible!')
        end
        col0 = col0 + kmis(end);
        pos0(2) = col0;
        
        %% Skip if the piece was already considered
        
        if sco >= -1
            k = find(all(POS == pos0,2),1);
            if isempty(k)
                POS = [POS; pos0];
                SCO = [SCO;  sco];
            else
                if SCO(k) < sco
                    return
                else
                    SCO(k) = sco;
                end
            end
        end
    end
end

%% Find possible scores

% Find deletions; col0 - end of the previous piece
N1 = sum(MStep == 2*Step + 1);
N2 = sum(MStep == 2*Step + 2);
sco = sco + N1 + N2;

% Do not allow insertions with length higher than piece length
if reverse
    ii = Row{1+N1} > row0 + col0 - ColR{1+N1};
else
    ii = Row{1+N1} < row0 - col0 + ColR{1+N1};
end

% It ALLOWS stay in the same column:
% ii = find(ColL{1+N1} <= col0+1 & ColR{1+N1} >= col0 & Row{1+N1} ~= row0);

% It DOES NOT ALLOWS stay in the same column:
ii = find(ii & ColL{1+N1} <= col0+1 & ColR{1+N1} > col0 & ...
    Row{1+N1} ~= row0);

if isempty(ii)
    return
end
% Score
dRow = Row{1+N1}(ii)-row0;
LLii = LL{1+N1}(ii);
L0 = repmat(L0Del,size(dRow));
L0(dRow>0) = L0Ins;
S = sco + abs(dRow) .* (EpsTol + exp(L0-LLii));
% S = sco + sinh(abs(dRow)./LLii) .* exp(1-LLii./L0);
% Remove distant pieces
jj = S < 0;
S = S(jj);
if isempty(S)
    return
end
ii = ii(jj);

%% Sort scores and save right sides for future modification

[S, jj] = sort(S);
ii = ii(jj);
Rowii = Row{1+N1}(ii);
ColRii = ColR{1+N1}(ii);

%% Shift right side of each piece

if N2 ~= 0
    for i = 1:length(S)
        kmis = find(~M(Rowii(i),ColRii(i)+2:end),N2,'first');
        if length(kmis) ~= N2 % End is found
            row = Rowii(i);
            colL = col0 + 1;
            colR = Ns;
            S = S(i);
            return
        end
        ColRii(i) = ColRii(i) + kmis(end);
    end
elseif any(ColRii==Ns)
    % End is found
    i = find(ColRii==Ns,1,'first');
    row = Rowii(i);
    colL = col0 + 1;
    colR = Ns;
    S = S(i);
    return
end

%% Remove pieces which were already considered

if sco >= -1
    pos = [Rowii, ColRii];
    [~,iPOS,ipos] = intersect(POS,pos,'rows');
    inew = setdiff((1:length(S)).',ipos);
    if ~isempty(inew)
        POS = [POS; pos(inew,:)];
        SCO = [SCO; S(inew)];
    end
    if ~isempty(iPOS)
        ii = S(ipos) < SCO(iPOS);
        if any(ii)
            SCO(iPOS(ii)) = S(ipos(ii));
        end
        ii = ~ii;
        if any(ii)
            jj = ipos(ii);
            S(jj) = [];
            if isempty(S)
                return
            end
            Rowii(jj) = [];
            ColRii(jj) = [];
        end
    end
end

%% Loop

Stepi = zeros(size(S));
for i = 1:length(S)
    [row,colL,colR,Stepi(i),POS,SCO,St] = SinglePiece(L0Del,L0Ins,S(i), ...
        Rowii(i),ColRii(i),Row,ColL,ColR,LL,Ns,POS,SCO,Step,MStep,M, ...
        reverse,EpsTol);
    if ~isempty(row)
        row = [Rowii(i); row]; %#ok
        colL = [col0+1; colL]; %#ok
        colR = [ColRii(i); colR]; %#ok
        S = St;
        return
    end
end
Step = max(Step,max(Stepi));