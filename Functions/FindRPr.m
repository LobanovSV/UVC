function Pr = FindRPr(Het,A,ALR,SEQ,QUAL,POSS,PAIR,LEFT,maxSL,dPos, ...
    dPOS,Qmax)

%% Parameters

Fact = -0.5 * log10(exp(1)) / maxSL;
NEQtol = 3;
maxErr = NEQtol * Qmax;

%% Find left and right alleles

NA = length(A);
Np = max(PAIR);
Err = zeros(Np,NA,Het);
POS = zeros(Np,4,NA,Het);
MinL = inf;
for j = 1:NA
    for i = 1:Het
        allele = [ALR{i,1},A{j},ALR{i,2}];
        [Err(:,j,i),POS(:,:,j,i)] = SimplePairAlignment( ...
            allele,SEQ,QUAL,POSS,PAIR,LEFT,[],[],maxSL,dPos,dPOS,[],Qmax);
    end
    MinL = min(MinL,length(A{j}));
end
In = reshape(any(Err ~= Err(:,1,:),2),Np,Het);

%% Replace maximum error for each read by SeqErr

MinErr = min(Err,[],[2,3]);
Err = Err - MinErr;
W = exp(-MinErr/maxErr);

%% Find positions for depth difference

PP = cell(Het,NA);
for j = 1:NA
    pp = round(linspace(1,length(A{j}),MinL));
    for i = 1:Het
        p2 = length(ALR{i,1});
        p1 = p2 - maxSL;
        p3 = p2 + length(A{j}) + 1;
        p4 = p3 + maxSL;
        PP{i,j} = [p1:p2,p2+pp,p3:p4];
    end
end

%% Find Pr

Pr = zeros(NA);
for i = 1:NA
    for j = 1:NA
        if Het == 1 && j < i
            Pr(i,j) = Pr(j,i);
            continue
        end
        E1 = Err(:,i,1);
        E2 = Err(:,j,Het);
        
        % First depth coverage
        ii = find(In(:,1) & E1<E2);
        [~,jj] = unique(POS(ii,[1,4],i,1),'rows');
        pos = POS(ii(jj),:,i,1);
        pp = PP{1,i};
        dc1 = zeros(size(pp));
        for k = 1:size(pos,1)
            for l = [1,3]
                if pos(k,l)<=pp(end) && pos(k,l+1)>=pp(1)
                    ii = pp>=pos(k,l) & pp<=pos(k,l+1);
                    dc1(ii) = dc1(ii) + 1;
                end
            end
        end
        
        % Second depth coverage
        ii = find(In(:,Het) & E2<E1);
        [~,jj] = unique(POS(ii,[1,4],j,Het),'rows');
        pos = POS(ii(jj),:,j,Het);
        pp = PP{Het,j};
        dc2 = zeros(size(pp));
        for k = 1:size(pos,1)
            for l = [1,3]
                if pos(k,l)<=pp(end) && pos(k,l+1)>=pp(1)
                    ii = pp>=pos(k,l) & pp<=pos(k,l+1);
                    dc2(ii) = dc2(ii) + 1;
                end
            end
        end
        
        Pr(i,j) = Fact * nansum((dc1 - dc2).^2 ./ (dc1 + dc2)) - ...
            sum(W.*min(maxErr,min(E1,E2)))/10;
    end
end

%% Convert to probabilities

Pr = Pr - max(Pr,[],'all');
Pr = 10 .^ Pr;

%% Normalize so that sum(sum(Pr)) = 1

% Pr = Pr / sum(Pr,'all');