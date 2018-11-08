function pValue = multipleHotelling(bin)

pValue=nan(length(bin.x),length(bin.y));

p1=hotelling(bin,bin,1,4);
p2=hotelling(bin,bin,2,5);
p3=hotelling(bin,bin,3,6);

p4=hotelling(bin,bin,1,2);
p5=hotelling(bin,bin,1,3);
p6=hotelling(bin,bin,1,5);
p7=hotelling(bin,bin,1,6);

p8=hotelling(bin,bin,2,3);
p9=hotelling(bin,bin,2,4);
p10=hotelling(bin,bin,2,6);

p11=hotelling(bin,bin,3,4);
p12=hotelling(bin,bin,3,5);

p13=hotelling(bin,bin,4,5);
p14=hotelling(bin,bin,4,6);

p15=hotelling(bin,bin,5,6);

for i=1:length(bin.x)
    for j=1:length(bin.y)

        pValue(i,j)=min([p1(i,j),p2(i,j),p3(i,j),p4(i,j),p5(i,j),...
            p6(i,j),p7(i,j),p8(i,j),p9(i,j),p10(i,j),p11(i,j),p12(i,j),p13(i,j),p14(i,j),p15(i,j)]);
        
    end
end

end

