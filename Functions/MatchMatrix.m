function [M0,pos0] = MatchMatrix(REF,Pos,dPOS,maxSL)

%% Shift left position

NREF = length(REF);

%% Create match matrix

i1 = max(1,min(NREF,Pos(1)) - dPOS);
i2 = min(NREF,Pos(end) + dPOS + maxSL - 1);
pos0 = i2 - maxSL + 1;
N = pos0 + 1 - i1;

REF = REF(i2:-1:i1).';

M0 = repmat(blanks(maxSL),N,1);
for i = 1:maxSL
    M0(:,i) = REF(end+2-i-N:end+1-i);
end