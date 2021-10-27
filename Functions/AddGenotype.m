function Genotype = AddGenotype(G,Gp,isBad)

N = size(G,1);
if nargin < 3
    isBad = false(N,1);
end

if isempty(G)
    NR = 0;
else
    NR = size(G,3);
end

Genotype = cell(N,1);
for i = 1:N
    genotype = cell(1,NR);
    for j = 1:NR
        genotype{j} = [int2str(G(i,1,j)) '/' int2str(G(i,2,j)) ...
            '(' NUM2STR(Gp(i,j)) '\%)'];
    end
    if isBad(i)
        genotype{1,end+1} = 'might be poorly aligned'; %#ok
    end
    if isempty(genotype)
        Genotype{i} = '';
    else
        Genotype{i} = {strjoin(genotype),', '};
    end
end
end

function y = NUM2STR(x)
y = num2str(100*x,'%.1f');
if y(end) == '0' && y(end-1) == '.'
    y = y(1:end-2);
end
end