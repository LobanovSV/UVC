function [G,Gp,Ap] = FindGp(G,A,RS,RE,RF,SEQ,QUAL,POSS,POSE,PAIR,LEFT, ...
    dPos,dPOS,RAM,FILE,i)

%% Parameters

NStep = 10;
ReCalcDistTol = 100;

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

Calc = true(1,NR);
Gp = zeros(1,NR);
ALR = cell(2,2);
GLoop = zeros(NStep,2,NR);
Ap = cell(NR,1);
for Step = 1:NStep
    if any(all(GLoop(1:Step-1,:,:)==G,[2,3]))
        break
    end
    GLoop(Step,:,:) = G;
    for i = find(Calc)
        for j = 1:2
            ALR{j,1} = strjoin(RF(j,1:2*i-1),'');
            ALR{j,2} = strjoin(RF(j,2*i+1:end),'');
        end
        Het = 1 + any(diff(G(1,:,[1:i-1,i+1:end]),1,2)~=0);
        Pr = FindRPr(Het,A{i},ALR,SEQ,QUAL,POSS,PAIR,LEFT,maxSL,dPos, ...
            dPOS,Qmax);
        
        [Pr0,Ind] = max(Pr,[],'all','linear');
        i1 = G(1,1,i);
        i2 = G(1,2,i);
        if Pr0 > Pr(i1,i2)
            i1 = 1;
            i2 = i1;
            if Pr0 > Pr(i1,i2)
                NA = length(A{i});
                [i1,i2] = ind2sub([NA,NA],Ind);
            end
        end
        
        Pr = max(triu(Pr),tril(Pr).');
        Pr = Pr / sum(Pr,'all');
        Gp(i) = abs(diff(maxk(Pr(:),2)));
        
        if i1 ~= G(1,1,i) || i2 ~= G(1,2,i)
            if G(1,1,i) ~= i1
                G(1,1,i) = i1;
                RF{1,2*i} = A{i}{i1};
            end
            if G(1,2,i) ~= i2
                G(1,2,i) = i2;
                RF{2,2*i} = A{i}{i2};
            end
            Calc(RE>=RS(i)-ReCalcDistTol & RS<=RE(i)+ReCalcDistTol) = true;
        end
        
        % Find allele probabilities
        Ap1 = max(Pr,[],1).';
        Ap2 = max(Pr,[],2);
        Ap{i} = max(Ap1/sum(Ap1),Ap2/sum(Ap2));
        
        Calc(i) = false;
    end
    if ~any(Calc)
        break
    end
end