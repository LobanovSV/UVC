function [Allele,FPOS] = LoadAlleleFPOS(AlleleRAM,FPOSRAM,RAM,FILE,i)

if RAM
    Allele = AlleleRAM{i};
    FPOS = FPOSRAM{i};
else
    FileLoad = fullfile(DirPCorUNIX,'BAMs','Temp',FILE, ...
        [int2str(i) '.mat']);
    load(FileLoad,'Allele','FPOS')
end