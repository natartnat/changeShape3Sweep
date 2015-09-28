function [D]=Get_file()
[filename,pathname] = uigetfile({'*.png';'*.jpg';'*.tif';'*.*'}, 'Get Image file');
        fich1 = fullfile(pathname,filename);
        D=imread(fich1);
       
           
        
