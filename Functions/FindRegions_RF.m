function RF = FindRegions_RF(RS,RE,NR,FASTA,posW)

if NR == 0
    RF = {FASTA};
    return
end

RS_F = posW(1) - 2 + RS;
RE_F = posW(1) - 2 + RE;
RF = cell(1,2*NR+1);
RF{1} = FASTA(1:RS_F(1));
for i = 1:NR-1
    RF{2*i} = FASTA(RS_F(i)+1:RE_F(i)-1);
    RF{2*i+1} = FASTA(RE_F(i):RS_F(i+1));
end
RF{end-1} = FASTA(RS_F(end)+1:RE_F(end)-1);
RF{end} = FASTA(RE_F(end):end);