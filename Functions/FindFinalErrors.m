function [dEp,iPCRp,POSp] = FindFinalErrors(PCR,STR,SEQ,QUAL,POSS,PAIR, ...
    LEFT,Np,maxSL,dPos,dPOS)

Ep = inf(Np,2);
iPCRp = ones(Np,2);
POSp = repmat(inf(Np,1) .* [-1,-1,1,1],1,1,2);
RIGHT = ~LEFT;

for i = 1:2
    if length(PCR) < i
        break
    end
    Ni = length(PCR(i).seq);
    
    [Ei,POSSi,POSEi] = SimpleAlignment(PCR(i).seq{1},SEQ,QUAL,POSS, ...
        maxSL,dPos,dPOS);
    
    if Ni ~= 1
        % Repmat for Ni-1 PCRs
        Ei = repmat(Ei,1,Ni);
        POSSi = repmat(POSSi,1,Ni);
        POSEi = repmat(POSEi,1,Ni);
        
        % Shift positions
        for j = 1:length(PCR(i).pos)
            ii = POSSi(:,1) > PCR(i).pos(j);
            POSSi(ii,:) = POSSi(ii,:) + PCR(i).dN(j,:);
            POSEi(ii,:) = POSEi(ii,:) + PCR(i).dN(j,:);
        end
        
        for j = 2:Ni
            k = PCR(i).Ind(j);
            
            % Choose intersecting reads
            pos1 = STR(i).PosS(k) + abs(PCR(i).dN(k,j));
            pos2 = STR(i).PosE(k) - abs(PCR(i).dN(k,j));
            II = (POSSi(:,1) < pos1 & POSEi(:,1) > STR(i).PosE(k)) | ...
                (POSSi(:,1) < STR(i).PosS(k) & POSEi(:,1) > pos2);
            
            if any(II)
                [Ei(II,j),POSSi(II,j),POSEi(II,j)] = SimpleAlignment( ...
                    PCR(i).seq{j},SEQ(II),QUAL(II),POSS(II),maxSL, ...
                    dPos,dPOS);
            end
        end
    end
    
    Epi = zeros(Np,Ni);
    POSpi = repmat(POSp(:,:,i),1,1,Ni);
    
    ii = PAIR(LEFT);
    Epi(ii,:) = Ei(LEFT,:);
    POSpi(ii,1,:) = POSSi(LEFT,:);
    POSpi(ii,2,:) = POSEi(LEFT,:);
    ii = PAIR(RIGHT);
    Epi(ii,:) = Epi(ii,:) + Ei(RIGHT,:);
    POSpi(ii,3,:) = POSSi(RIGHT,:);
    POSpi(ii,4,:) = POSEi(RIGHT,:);
    
    Epi = Epi + PCR(i).dE.';
    
    [Ep(:,i),iPCRp(:,i)] = min(Epi,[],2);
    POSp(:,:,i) = POSpi(:,:,1);
    for j = 2:Ni
        ii = iPCRp(:,i) == j;
        if any(ii)
            POSp(ii,:,i) = POSpi(ii,:,j);
        end
    end
end

dEp = diff(Ep,1,2);