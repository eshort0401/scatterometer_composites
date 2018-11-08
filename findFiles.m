function [filePathArray] = findFiles(searchTerm)

p=path;

k=strfind(p,':');
pathList=cell(1,length(k));

pathList{1}=p(1:k(1)-1);
filePathArray={};
fileCount=1;

for i=2:length(pathList)
    pathList{i}=p(k(i-1)+1:k(i)-1);
end

% Find files.
for i=1:length(pathList)
    
    filePath=dir([pathList{i} '/' searchTerm]);
    filePath={filePath.name};
    
    if ~isempty(filePath)
        for j=1:length(filePath)
            filePathj=[pathList{i} '/' filePath{j}];
            filePathArray{fileCount}=filePathj;
            fileCount=fileCount+1;
        end
    end
        
end

end

