function [E,POSS,POSE] = SimpleAlignment(REF,SEQ,QUAL,POSS0,maxSL, ...
    dPos,dPOS,PAIR,LEFT)

if nargin < 8
    PAIR = [];
end
if nargin < 9
    LEFT = [];
end

%% Parameters

MinSL = 30;

%% Find the best alignment

[~,iiPos] = sort(POSS0);

N = length(SEQ);
E = zeros(N,1);
POSS = zeros(N,1);
POSE = zeros(N,1);
Pos = [0, 0];
for i = 1:length(iiPos)
    j = iiPos(i);
    poss = POSS0(j);
    if poss > Pos(2)
        Pos = poss + [0, dPos];
        [M0,pos0] = MatchMatrix(REF,Pos,dPOS,maxSL);
    end
    [E(j),POSS(j),POSE(j)] = FindBestPos(SEQ{j},QUAL{j},M0,pos0);
end

%% Correct pairs

if ~isempty(PAIR)
    % Find indexes of pairs
    Ind = zeros(max(PAIR),2);
    Indt = 1:length(PAIR);
    Ind(PAIR(LEFT),1) = Indt(LEFT);
    Ind(PAIR(~LEFT),2) = Indt(~LEFT);
    Ind(~all(Ind,2),:) = [];
    
    for i = 1:size(Ind,1)
        jl = Ind(i,1);
        jr = Ind(i,2);
        
        % Modify positions
        if POSE(jr) < POSS(jl) + MinSL
            [E(jl),POSS(jl),POSE(jl),E(jr),POSS(jr),POSE(jr)] = ...
                FindBestPairPos(POSS0(jl),SEQ{jl}, ...
                QUAL{jl},POSS0(jr),SEQ{jr},QUAL{jr},REF,dPOS,MinSL);
        end
        
        % Remove left nucleotides
        n = POSS(jl) - POSS(jr);
        if n > 0 && n < MinSL
            E(jr) = E(jr) - QUAL{jr}(1:n) * ...
                (SEQ{jr}(1:n) ~= REF(POSS(jr):POSS(jr)-1+n)).';
            if E(jr) < 0
                error('Wrong')
            end
        end
        
        % Remove right nucleotides
        n = POSE(jl) - POSE(jr);
        if n > 0 && n < MinSL
            E(jl) = E(jl) - QUAL{jl}(end+1-n:end) * ...
                (SEQ{jl}(end+1-n:end) ~= REF(POSE(jl)+1-n:POSE(jl))).';
            if E(jl) < 0
                error('Wrong')
            end
        end
    end
end