function [NPCR,Qmaxi] = FindPCRerr(NR,Qmax,ml)

%% Parameters

NPCRL = 4;

%% Calculate

NPCR = round(NR / NPCRL);

if nargin > 2
    Qmaxi = Qmax / NPCR / (1 + double(ml==1));
end