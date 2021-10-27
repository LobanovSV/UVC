% This script performs in-depth read alignment and variant calling.
% It calls any type of variant including SNP, INDEL, STR, QTR, and any
% their combination.
%
% Please cite the following publication if you use this code:
% Lobanov et. al., bioRxiv 2021.07.16.452643.
%
% Please read instructions at https://github.com/LobanovSV/UVC.git
% before using this script.
%
% Software is created by S. Lobanov in 2021.

clc
clear
SetPath

%% Parameters

% Choose folder
% It must contain the sub-folder SAMs with sam-files to be analysed
DirWork = uigetdir;

%% Default parameters

% Allele frequency tolerance (alleles with frequency below this value are
% considered artefacts and are removed)
AFTol = 1e-3;

% Save final alignment?
% Note! You need pdftex installed in order to save final alignment.
% If you don't have it installed, you should choose 'false'
SavePdf = true;

% Version of the reference genome which was used for read alignment
% The hg19 and hg38 genomes are provided and located in the 'Data'
% sub-folder. If you wish to use your own genome, look into structure of
% the hg19 and hg38 genomes and create your own in the same way
RefGenome = 'hg19'; % hg38

% Extend regions
dPos = 100;
dPOS = 200;

%% DirSave

DirLoad = fullfile(DirWork,'SAMs');
FileLoad = fullfile(DirLoad,'ID.mat');
load(FileLoad,'ID')

%% Load BAMs

[pos0,chr,posW,FASTA,ID,~,SEQ,QUAL,POSS,POSE,PAIR,LEFT] = LoadSAMs(...
    DirLoad,ID,dPOS,[],[],[],RefGenome);

%% Single Sample Alignment

N = length(ID);
Allele = cell(N,1);
Alg = cell(N,1);
Err = cell(N,1);
for i = 1:N
    clc
    fprintf('Single sample alignment: sample %u from %u\n',i,N)
    [Allele{i},Alg{i},Err{i},IIREM] = SingleSampleAlignment(FASTA, ...
        SEQ{i},QUAL{i},POSS{i},POSE{i},PAIR{i},LEFT{i},posW,dPos,dPOS);
    
    % Remove reads with bad score
    if ~isempty(IIREM)
        SEQ{i}(IIREM) = [];
        QUAL{i}(IIREM) = [];
        POSS{i}(IIREM) = [];
        POSE{i}(IIREM) = [];
        PAIR{i}(IIREM) = [];
        LEFT{i}(IIREM) = [];
        ID{i} = [ID{i} ', ' int2str(sum(IIREM)) ...
            ' reads were removed due to large fraction of mismatches'];
    end
end
Err = vertcat(Err{:});

%% Checks and removes

IIREM = cellfun(@length,Allele) ~= 4;
ID(IIREM) = [];
SEQ(IIREM) = [];
QUAL(IIREM) = [];
POSS(IIREM) = [];
POSE(IIREM) = [];
PAIR(IIREM) = [];
LEFT(IIREM) = [];
Allele(IIREM) = [];
Alg(IIREM) = [];
Err(IIREM) = [];

%% Joint Sample Alignment

[RS,RE,A,AisN,G,Gp,IDL,IDLp,Allele,Alg,Err] = JointSampleAlignment( ...
    Allele,Alg,Err,SEQ,QUAL,POSS,POSE,PAIR,LEFT,FASTA,posW,dPos,dPOS, ...
    AFTol);

%% Add genotype

Genotype = AddGenotype(G,Gp);

%% Save

clc
disp('Save results...')
FileSave = fullfile(DirWork,'Results.mat');
save(FileSave,'chr','ID','Allele','Alg','Err','Genotype', ...
    'RS','RE','A','AisN','G','Gp','IDL','IDLp','pos0','posW','AFTol', ...
    'IIREM','-v7.3')

%% Draw

if SavePdf
    clc
    disp('Create pdf...')
    SaveAlignment(ID,Allele,Alg,DirWork,'Alignment',Genotype);
end