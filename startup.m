% Add if statement to automatically detect which machine I'm on. 

% Personal Macbook.
if ismac
    addpath(genpath('/Volumes/Ewan''s Hard Drive/Data'));
    addpath(genpath('/Volumes/EWAN_HD_2/Data'));
    addpath(genpath('/Volumes/Ewan''s Hard Drive/Ewan_Short_Masters/MatLab Scripts'));
    savepath('/Volumes/Ewan''s Hard Drive/Ewan_Short_Masters/MatLab Scripts/Data Analysis/pathdef.m');
end
   
% Uni Unix box machines.
if isunix && not(ismac)
    username=char(java.lang.System.getProperty('user.name'));
    addpath(genpath(['/media/' ...
       username '/Ewan''s Hard Drive/Data']));
    addpath(genpath(['/media/' ...
       username '/EWAN_HD_2/Data']));
    addpath(genpath(['/media/' ...
        username '/Ewan''s Hard Drive/Ewan_Short_Masters/MatLab Scripts']));
    savepath(['/media/' ...
        username '/Ewan''s Hard Drive/Ewan_Short_Masters/MatLab Scripts/Data Analysis/pathdef.m']);
    clear username; 
end