%% function [ new_data ] = twoDBaseProfile(temp)
function [ new_data, range ] = twoDBaseProfile(temp)
figure('Name','detection','NumberTitle','off'),imshow(temp);
close detection;
[out ,range] = myCrop( temp );
line = ceil((range(4)-range(3))/2); %% down-up |y1-y2 %height
halfHeight = range(4)-line;
real_data = temp(halfHeight:range(4),range(1):range(2)); %% half lower
[real_data, range] = myCrop( real_data ); %% get bound
interpolation_data = fliplr(flipud(real_data)); %% flip half lower
%new_data = ([interpolation_data; real_data]);
%figure,imshow(new_data);
new_data = imfill([interpolation_data; real_data]);%% concat data

%% ===== mod
%new_data = imresize(new_data,[size(new_data,2) size(new_data,2)]); % symmetric
[new_data, range] = myCrop( new_data ); %% get bound
figure('Name','basedProfile','NumberTitle','off'),imshow(new_data);
close basedProfile;
end

