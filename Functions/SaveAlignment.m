function SaveAlignment(Names,Allele,Alg,Dir,Name,AddText,STR,QTR,Pos, ...
    Col,ColD)
if nargin < 6 || isempty(AddText)
    AddText = cell(size(Alg));
end
if nargin < 7
    STR = [];
end
if nargin < 8
    QTR = [];
end
if nargin < 9 || isempty(Pos)
    Pos = cell(size(Alg));
end
if nargin < 10
    Col = [];
end
if nargin < 11
    ColD = [];
end

%% Parameters

LetterWidth = 8.7;
ColWidth = 0.01;
ColDWidth = 0.1;
LetterHeight = 12;
NMaxLines = 900;
TextType = '*ACGT';
Style = {'\color{black}','\color{green}','\color{blue}', ...
    '\color{BurntOrange}','\color{red}','\color{black}'};
DoNotMark = ' <>';

%% Remove empty alleles

ii = cellfun(@isempty,Allele);
if all(ii)
    return
elseif any(ii)
    Names(ii) = [];
    Allele(ii) = [];
    Alg(ii) = [];
    AddText(ii) = [];
    Pos(ii) = [];
end

%% Remove some lines if there are too many lines

for i = 1:length(Alg)
    NL = cellfun(@(x) size(x,1),Alg{i});
    NLrem = max(0,NL - NMaxLines / length(NL));
    NLrem = max(0,round((sum(NL) - NMaxLines) * NLrem / sum(NLrem)));
    if ~isnan(NLrem(1)) && any(NLrem ~= 0)
        for j = 1:length(NLrem)
            if NLrem(j) ~= 0
                ii = randperm(NL(j),NLrem(j));
                Alg{i}{j}(ii,:) = [];
            end
        end
        Names{i} = sprintf('%s, %u reads are not shown',Names{i}, ...
            sum(NLrem));
    end
end

%% Additional parameters

N = length(Allele);

% Page size
paperwidth = 0;
paperheight = 0;
for i = 1:N
    ph = max(cellfun(@(x) size(x,2),Alg{i}));
    paperwidth = max(paperwidth,ph);
    ph = size(AddText{i},1) + ...
        sum(cellfun(@(x) size(x,1),Allele{i})) + ...
        sum(cellfun(@(x) size(x,1),Alg{i}));
    paperheight = max(paperheight,ph);
end
paperwidth = paperwidth * LetterWidth + sum(Col) * ColWidth + ...
    sum(ColD) * ColDWidth;
paperwidth = int2str(paperwidth);
paperheight = int2str((paperheight+2)*LetterHeight);

clear ph

%% Create temporary folder for TeX

FileTeX = fullfile(Dir,[Name '.tex']);

TeX = fopen(FileTeX,'W');

%% Preambula

fprintf(TeX,'%s\n\n','\documentclass[10pt]{article}');

fprintf(TeX,'%s\n','\usepackage[table,dvipsnames]{xcolor}');
fprintf(TeX,'%s%s%s%s%s%s\n','\usepackage[top=0pt,bottom=0pt,', ...
    'left=0pt,right=0pt,paperwidth=',paperwidth,'pt,paperheight=', ...
    paperheight,'pt]{geometry}');
fprintf(TeX,'%s\n','\setlength\parindent{0pt}');
fprintf(TeX,'%s%s%s\n\n','\setlength{\tabcolsep}{0pt}');

fprintf(TeX,'%s\n\n','\begin{document}');


%% Single alignment

for i = 1:N
    SingleAlignment(Names{i},Allele{i},Alg{i},TeX,LetterWidth, ...
        LetterHeight,TextType,Style,DoNotMark,AddText{i},STR,QTR, ...
        Pos{i},Col,ColD,ColWidth,ColDWidth)
end

%% Finalize TeX and run pdfLaTeX

fprintf(TeX,'%s\n','\end{document}');
fclose(TeX);

%% Remove temporary folder, close files, etc.

CD = cd;
cd(Dir)
% system('pdflatex Tex.tex');
if ispc
    Options = '--extra-mem-bot=1000000000 --extra-mem-top=1000000000 ';
else
    Options = '';
end
system(['pdflatex --enable-write18 ' Options Name '.tex']);
cd(CD)
FileDel = fullfile(Dir,[Name '.tex']);
delete(FileDel)
FileDel = fullfile(Dir,[Name '.log']);
delete(FileDel)
FileDel = fullfile(Dir,[Name '.aux']);
delete(FileDel)
% movefile(fullfile(FileSave,'Tex.pdf'),'Tex.pdf')
% rmdir(FileSave,'s')
% open('Tex.pdf')
end

function SingleAlignment(Name,Allele,Alg,TeX,LetterWidth,LetterHeight, ...
    TextType,Style,DoNotMark,AddText,STR,QTR,Pos,Col,ColD,ColWidth, ...
    ColDWidth)

N2 = max(cellfun(@(x) size(x,2),Alg));

%% Replace underlines

ii = find(Name=='_');
for i = length(ii):-1:1
    Name = [Name(1:ii(i)-1) '\' Name(ii(i):end)];
end

%% Page size

paperwidth = N2 * LetterWidth + sum(Col) * ColWidth + ...
    sum(ColD) * ColDWidth;
paperwidth = int2str(paperwidth);
ph = sum(cellfun(@(x) size(x,1),Allele));
paperheight = int2str((2*ph + size(AddText,1) + ...
    sum(cellfun(@(x) size(x,1),Alg)))*LetterHeight);

clear ph

%% Find Positions

DrawCol = ~isempty(Col) || ~isempty(ColD) || ...
    (~isempty(QTR) && ~isempty(QTR.Motif)) || ...
    (~isempty(STR) && ~isempty(STR.Motif));
if DrawCol && isempty(Pos)
    ii = Allele{end} ~= '*';
    Pos = NaN(size(Allele{end}));
    Pos(ii) = 1:sum(ii);
end

%% Write text

fprintf(TeX,'%s%s%s%s%s\n\n','\eject \pdfpagewidth=',paperwidth, ...
    'pt \pdfpageheight=',paperheight,'pt');

fprintf(TeX,'%s%s%s\n\n','\textbf{',Name,'}');

for i = 1:size(AddText,1)
    if i == 1 || ~strcmp(AddText{1},AddText{i})
        fprintf(TeX,'%s%s%s\n\n','\textbf{',AddText{i},'}');
    end
end

%% Find columns and write them to the file

Columns = repmat('c',1,N2);
if ~isempty(QTR) && ~isempty(QTR.Motif)
    rl = length(QTR.Motif{1}{1});
    Col = unique([Col,QTR.PosS-1:rl:QTR.PosE]);
end
if ~isempty(STR) && ~isempty(STR.Motif)
    rl = length(STR.Motif{1});
    ColD = unique([ColD,STR.PosS-1,STR.PosE-mod(STR.PosE+1-STR.PosS,rl)]);
end
if DrawCol
    % Single columns
    C = zeros(2,length(Col));
    for i = 1:length(Col)
        C(1,i) = find(Pos==Col(i),1);
        C(2,i) = find(Pos==Col(i)+1,1) - 1;
    end
    C = unique(C);
    % Double columns
    dC = zeros(length(ColD),1);
    for i = 1:length(ColD)
        dC(i) = find(Pos==ColD(i),1);
    end
    % Merge
    C = unique([C(:);dC(:)]);
    Double = ismember(C,dC);
    
    for i = length(C):-1:1
        j = C(i);
        if Double(i)
            Vert = '||';
        else
            Vert = '|';
        end
        Columns = [Columns(1:j) Vert Columns(j+1:end)];
    end
end
fprintf(TeX,'%s%s%s%s%s\n','\begin{tabular}{',Columns,'}');

%% Write alleles and alignment

for i = 1:length(Allele)
    if ~isempty(Allele{i})
        Line = AddColor(Allele{i},Allele{i},Allele{end},TextType, ...
            Style,DoNotMark,true);
        fprintf(TeX,'%s\n',Line);
    end
    if ~isempty(Alg{i})
        Line = AddColor(Alg{i},Allele{i},Allele{end},TextType, ...
            Style,DoNotMark,false);
        if iscell(Line)
            fprintf(TeX,'%s\n',strjoin(Line,'\n'));
        else
            fprintf(TeX,'%s\n',Line);
        end
    end
end

% fprintf(TeX,'%s%s%s\n',repmat('\color{white}G&',1,N2-1), ...
%     '\color{white}G\\');

% fprintf(TeX,'%s\n',' & & &\cellcolor{red}  & &C&T&C&A');
fprintf(TeX,'%s\n\n','\end{tabular}');
end

function AlgOut = AddColor(AlgChar,LocAllele,RefAllele,TextType, ...
    Style,DoNotMark,IsAlele)
N1 = size(AlgChar,1);
Alg = num2cell(AlgChar);

% Replace special symbols
for i = 1:numel(Alg)
    if strcmp(Alg{i},'<')
        Alg{i} = '\textless';
    end
    if strcmp(Alg{i},'>')
        Alg{i} = '\textgreater';
    end
end

NonEmpty = true(size(AlgChar));
for i = 1:length(DoNotMark)
    NonEmpty = NonEmpty & AlgChar ~= DoNotMark(i);
end
% Gray
ii = NonEmpty & AlgChar == RefAllele;
if IsAlele
    AddText = '\cellcolor{gray}\color{white}\bf{}';
else
    AddText = '\color{gray}';
end
Alg(ii) = arrayfun(@(x) {[AddText x]}, AlgChar(ii));
NonEmpty = NonEmpty & ~ii;
% New nucleotides
for i = 1:length(TextType)
    ii = NonEmpty & AlgChar == TextType(i);
    if IsAlele
        AddText = ['\cell' Style{i}(2:end) '\color{white}\bf{}'];
    else
        AddText = Style{i};
    end
    Alg(ii) = arrayfun(@(x) {[AddText x]}, AlgChar(ii));
    NonEmpty = NonEmpty & ~ii;
end
% The rest nucleotides
if IsAlele
    AddText = ['\cell' Style{end}(2:end) '\color{white}\bf{}'];
else
    AddText = Style{end};
end
Alg(NonEmpty) = arrayfun(@(x) {[AddText x]}, AlgChar(NonEmpty));

if ~IsAlele
    ii = AlgChar ~= LocAllele;
    for i = 1:length(DoNotMark)
        ii = ii & AlgChar ~= DoNotMark(i);
    end
    Alg(ii) = arrayfun(@(x) {['\cellcolor{Lavender}' x]}, AlgChar(ii));
end

Alg(:,1:end-1) = cellfun(@(x) {[x ' & ']}, Alg(:,1:end-1));

Alg = [Alg, repmat({' \\ '},N1,1)];

if N1 == 1
    AlgOut = cell2mat(Alg);
else
    AlgOut = cell(N1,1);
    for i = 1:N1
        AlgOut{i} = cell2mat(Alg(i,:));
    end
end
end