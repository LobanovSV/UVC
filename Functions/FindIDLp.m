function [IDL,IDLp] = FindIDLp(G,A,IDL_L,IDL_Ind,RF,SEQ,QUAL,POSS,POSE, ...
    PAIR,LEFT,dPos,dPOS,RAM,FILE,i)

%% Parameters

%% Load data

if ~RAM
    File = fullfile(DirPCorUNIX,'BAMs','Temp',FILE,[int2str(i) '.mat']);
    load(File,'SEQ','QUAL','PAIR','LEFT','POSS','POSE')
end

%% Maximum quality

Qmax = max(cellfun(@max,QUAL));
maxSL = 1 + max(POSE - POSS);

%% Fill RF with G & A, find unique lengths

RF = repmat(RF,2,1);
NR = length(A);
for i = 1:NR
    for j = 1:2
        RF{j,2*i} = A{i}{G(1,j,i)};
    end
end

%% The main loop

IDL = zeros(8,NR);
IDLp = zeros(8,NR);
ALR = cell(2,2);
for i = 1:NR
    if isempty(IDL_L{1,i})
        continue
    end
    
    for j = 1:2
        ALR{j,1} = strjoin(RF(j,1:2*i-1),'');
        ALR{j,2} = strjoin(RF(j,2*i+1:end),'');
    end
    Het = 1 + any(diff(G(1,:,[1:i-1,i+1:end]),1,2)~=0);
    Pr = FindRPr(Het,A{i},ALR,SEQ,QUAL,POSS,PAIR,LEFT,maxSL,dPos, ...
        dPOS,Qmax);
    Pr = Pr(:);
    
    for j = 1:8
        if isempty(IDL_L{j,i})
            continue
        end
        n = length(IDL_L{j,i});
        pr = zeros(n,1);
        for k = 1:n
            pr(k) = max(Pr(IDL_Ind{j,i}(:,k)));
        end
        % Normalize so that sum(pr) = 1
        pr = pr / sum(pr);
        [pr2,i2] = maxk(pr,2);
        IDL(j,i) = IDL_L{j,i}(i2(1));
        IDLp(j,i) = pr2(1) - pr2(2);
    end
end