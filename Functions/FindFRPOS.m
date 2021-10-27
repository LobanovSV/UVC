function [FPOS,RPOS,RLL,NN2] = FindFRPOS(Allele,RAM,FILE,i)
% FPOS and RPOS are non-inclusive intervals (RS,RE)
% FPOS are (0,1,2,...,end-1,end,end+1)
% RPOS include left and right nucleotides i-1,i,i+1

%% Load data

if ~RAM
    File = fullfile(DirPCorUNIX,'BAMs','Temp',FILE,[int2str(i) '.mat']);
    load(File,'Allele')
end

%% FPOS

FPOS = [true,Allele{end}~='*',true];

%% RPOS

N2 = length(Allele{end});
rpos = zeros(1,N2+2);
rpos(FPOS) = 1:sum(FPOS);
FPOS = find(FPOS) - 1;

NN2 = 2 + sum(Allele{end}~='*');

if ~RAM
    save(File,'FPOS','-append')
    FPOS = [];
end

if nargout < 2
    return
end

%% RPOS

RPOS = cell(1,2);
RLL = cell(1,2);
for i = 1:2
    if ~isempty(Allele{i})
        Diffpr = zeros(1,0);
        Diff = find(Allele{i}~=Allele{end});
        Diff = unique([Diff,Diff+1,Diff+2]);
        while ~isequal(Diffpr,Diff)
            Diffpr = Diff;
            Diff0 = Diff(rpos(Diff)==0);
            Diff = unique([Diff,Diff0-1,Diff0+1]);
        end
        if isempty(Diff)
            RPOS{i} = zeros(0,1);
            RLL{i} = zeros(0,1);
            continue
        end
        ii = diff(Diff) ~= 1;
        iL = Diff([true,ii]);
        iR = Diff([ii,true]);
        RPOS{i} = rpos(Diff).';
        RPOS{i}(RPOS{i}==0) = [];
        RLL{i} = zeros(size(RPOS{i}));
        NotAst = [false,xor(Allele{i}=='*',Allele{end}=='*'),false];
        for j = 1:length(iL)
            LL = sum(NotAst(iL(j):iR(j)));
            RLL{i}(RPOS{i}>=rpos(iL(j)) & RPOS{i}<=rpos(iR(j))) = LL;
        end
    end
end