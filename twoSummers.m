% load('twoSummerFilenamesDatecells.mat')
% 
% bin=binASCAT(ASCATfiles,.25,-13,6,.25,94,130,16,dateCellASCAT);
% save('ASCATtwoSummerBor','bin','-v7.3');
% clear bin
% 
% bin=binASCAT(ASCATfiles,.25,-13,6,.25,130,152,16,dateCellASCAT);
% save('ASCATtwoSummerPap','bin','-v7.3');
% clear bin
% 
% bin=binHY2SCAT(HY2SCATfiles,.25,-13,6,.25,94,130,12,dateCellHY2SCAT);
% save('HY2SCATtwoSummerBor','bin','-v7.3');
% clear bin
% 
% bin=binHY2SCAT(HY2SCATfiles,.25,-13,6,.25,130,152,12,dateCellHY2SCAT);
% save('HY2SCATtwoSummerPap','bin','-v7.3');
% clear bin
 
% bin=binRSCATorOSCAT(OSCATfiles,.25,-13,6,.25,94,130,12,dateCellOSCAT,'OSCAT');
% save('OSCATtwoSummerBor','bin','-v7.3');
% clear bin

% bin=binRSCATorOSCAT(OSCATfiles,.25,-13,6,.25,130,152,12,dateCellOSCAT,'OSCAT');
% save('OSCATtwoSummerPap','bin','-v7.3');
% clear bin

bin=binWRF(WRFfiles,.25,-13,6,.25,94,130,6,dateCellWRF,'WRF');
save('WRFtwoSummerBor','bin','-v7.3');
clear bin

bin=binWRF(WRFfiles,.25,-13,6,.25,130,152,6,dateCellWRF,'WRF');
save('WRFtwoSummerPap','bin','-v7.3');
clear bin