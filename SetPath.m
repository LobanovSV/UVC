function SetPath

Path = {fullfile(cd,'Functions')};
for i = 1:length(Path)
    if ~contains(path,Path{i})
        path(Path{i},path)
    end
end