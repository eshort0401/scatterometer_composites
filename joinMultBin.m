function bin = joinMultBin(binCell)

bin=joinBin(binCell{1},binCell{2});

for i=3:length(binCell)
   
    bin=joinBin(bin,binCell{i});
    
end

fprintf('Joining bins. \n');

end

