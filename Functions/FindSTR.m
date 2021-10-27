function [Motif,ML,PosS,PosE] = FindSTR(Fasta,pos,MinNR)

%% Default parameters

if nargin < 2 || isempty(pos)
    pos = [1, length(Fasta)];
end
if nargin < 3 || isempty(MinNR)
    MinNR = 8;
end

%% Parameters

MaxML = 20;

%% Find all repeats

mm = MinNR - 1;
ML = [];
PosS = [];
PosE = [];
for ml = 1:MaxML
    Nref = pos(2) + 1 - pos(1) - ml;
    
    Match = find(Fasta(pos(1):pos(2)-ml) == Fasta(pos(1)+ml:pos(2)));
    if isempty(Match)
        continue
    end
    % Left and right boundaries
    LB = Match(diff([-1, Match]) ~= 1);
    RB = Match(diff([Match, Nref+2]) ~= 1);
    
    % Long matches
    ii = RB + 1 - LB >= mm * ml;
    if any(ii)
        LB = LB(ii).';
        RB = RB(ii).'+ml;
        % Check if STR is near boundary
%         if (any(LB==1)&&pos(1)~=1) || ...
%                 (any(RB==Nref)&&pos(2)~=length(Fasta))
%             error('Consider')
%         end
        RB = RB - mod(RB+1-LB,ml);
        
        PosS = [PosS; LB]; %#ok
        PosE = [PosE; RB]; %#ok
        ML = [ML; repmat(ml,length(LB),1)]; %#ok
    end
end

%% Check wether some repeat is found

if isempty(PosS)
    Motif = cell(0);
    return
end

%% Shift by pos(1)-1

PosS = PosS + pos(1)-1;
PosE = PosE + pos(1)-1;

%% Remove doublets, etc.

[PosS, ii] = unique(PosS);
PosE = PosE(ii);
ML = ML(ii);

%% Find motifs

NML = length(ML);
Motif = cell(NML,1);
for i = 1:NML
    Motif{i} = Fasta(PosS(i):PosS(i)-1+ML(i));
end

%% Remove non-A,C,G,T

iRem = false(NML,1);
for i = 1:NML
    iRem(i) = any(Motif{i}~='A' & Motif{i}~='C' & Motif{i}~='G' & ...
        Motif{i}~='T');
end
if any(iRem)
    Motif(iRem) = [];
    ML(iRem) = [];
    PosS(iRem) = [];
    PosE(iRem) = [];
end
% NML = length(ML);

%% Check wether there are intersections

% % Find intersections
% IntSct = (PosS <= PosS.' & PosS.' <= PosE) | ...
%     (PosS <= PosE.' & PosE.' <= PosE);
% IntSct(1:NML+1:end) = false;
% [i1, i2] = find(IntSct);
% clear IntSct
% 
% % Remove repetitions
% i12 = [i1, i2];
% i12 = unique(sort(i12,2),'rows');
% i1 = i12(:,1);
% i2 = i12(:,2);
% clear i12
% 
% if ~issorted(PosS)
%     error('PosS must be sorted here')
% end
% 
% % Replace
% for i = 1:length(i1)
%     NR1 = (PosE(i1(i)) + 1 -PosS(i1(i))) / ML(i1(i));
%     NR2 = (PosE(i2(i)) + 1 -PosS(i2(i))) / ML(i2(i));
%     if NR2 >= NR1 && PosE(i1(i)) < PosS(i2(i)) + ML(i2(i))
%         PosS(i2(i)) = PosS(i2(i)) + ML(i2(i));
%     elseif NR1 >= NR2 && PosE(i1(i)) - ML(i1(i)) < PosS(i2(i))
%         PosE(i1(i)) = PosE(i1(i)) - ML(i1(i));
%     else
%         disp(Fasta(min(PosS(i1(i)),PosS(i2(i))):max(PosE(i1(i)), ...
%             PosE(i2(i)))))
%         disp(Fasta(PosS(i1(i)):PosE(i1(i))))
%         disp(Fasta(PosS(i2(i)):PosE(i2(i))))
%         error('Consider')
%     end
% end
% 
% if any(sum(IntSct)~=1)
%     error('Consider')
% end