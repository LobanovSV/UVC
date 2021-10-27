function [A,AisN,RS,RE] = RemoveTheSameSeq(A,AisN,RS,RE)

%% Forward

NA = length(A);
Break = false;
for j = 1:length(A{1})
    if any(AisN(:,j)~=1)
        Break = true;
        break
    end
    c = A{1}(j);
    for k = 2:NA
        if length(A{k})<j || A{k}(j) ~= c
            Break = true;
            break
        end
    end
    if Break
        break
    end
end
j = j - Break;
if j ~= 0
    AisN(:,1:j) = [];
    for k = 1:NA
        A{k}(1:j) = [];
    end
    RS = RS + j;
end

%% Reverse

if any(AisN(:,end)~=0)
    return
end

Break = false;
for j = 1:length(A{1})
    if any(AisN(:,end+1-j)~=1)
        Break = true;
        break
    end
    c = A{1}(end+1-j);
    for k = 2:NA
        if length(A{k})<j || A{k}(end+1-j) ~= c
            Break = true;
            break
        end
    end
    if Break
        break
    end
end
j = j - Break;
if j ~= 0
    AisN(:,end-j:end-1) = [];
    for k = 1:NA
        A{k}(end+1-j:end) = [];
    end
    RE = RE - j;
end