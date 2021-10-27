function [IDL_L,IDL_Ind] = Find_IDL_Ind(A)
%  The order: Sum, Max, Min, Diff, Add3, Dom3, Rec3, Het3

NR = length(A);

IDL_L = cell(8,NR);
IDL_Ind = cell(8,NR);
for i = 1:NR
    l = cellfun(@length,A{i});
    l = l - l(1);
    if all(l==0)
        continue
    end
    l3 = mod(l,3);
    N = length(l)^2;
    
    for j = 1:8
        if j > 4 && all(l3==0)
            continue
        end
        switch j
            case 1
                ll = l + l.';
            case 2
                ll = max(l,l.');
            case 3
                ll = min(l,l.');
            case 4
                ll = abs(l - l.');
            case 5
                ll = l3 + l3.';
            case 6
                ll = max(l3,l3.');
            case 7
                ll = min(l3,l3.');
            case 8
                ll = abs(l3 - l3.');
        end
        [llu,~,ii] = unique(ll(:));
        Nu = length(llu);
        if Nu ~= 1
            IDL_L{j,i} = llu;
            IDL_Ind{j,i} = false(N,Nu);
            for k = 1:Nu
                IDL_Ind{j,i}(:,k) = ii == k;
            end
        end
    end
end