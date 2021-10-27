function [REF,posW] = AddInDel(REF,INDEL,posW)

if isempty(INDEL)
    return
end
ispos = nargin > 2 && nargout > 1;

%% Replace

N = size(INDEL,1);
for i = N:-1:1
    p1 = INDEL(i,1);
    p2 = INDEL(i,2);
    dp = p2 - p1 - 1;
    if p1 < p2
        REF = REF([1:p1,p2:end]);
        if ispos
            if posW(1) > p1
                posW(1) = max(p1,posW(1)-dp);
            end
            if posW(2) > p1
                posW(2) = max(p1+1,posW(2)-dp);
            end
        end
    else
        REF = REF([1:p1,p2:p1,p1+1:end]);
        if ispos
            if posW(1) > p1
                posW(1) = posW(1) - dp;
            elseif posW(1) > p2
                posW(1) = p2;
            end
            if posW(2) > p1
                posW(2) = posW(2) - dp;
            elseif posW(2) >= p2
                posW(2) = p1 - dp;
            end
        end
    end
end