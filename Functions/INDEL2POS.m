function [pS,pE] = INDEL2POS(INDEL)
% Find start and end of INDEL

%% Insertion

pS = INDEL(:,2);
pE = INDEL(:,1);

%% Deletion

ii = pE < pS;
pS(ii) = INDEL(ii,1) + 1;
pE(ii) = INDEL(ii,2) - 1;