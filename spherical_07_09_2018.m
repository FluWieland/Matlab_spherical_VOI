%% clear workspace
%also keep the run_count variable. this variable counts how many times you
%have run the script since you opened MATLAB last time, so no unnecessary
%questions are asked twice

clearvars -except counter coordinate_means_fiji;

%%
folder_name = inputdlg('Please enter the name of the folder to be created');

%%
message_folder = strcat('Please select the path to the folder containing the images. I',...
    'n this folder a subfolder will be created, containing the cropped images.');
uiwait(msgbox(message_folder))

%% get paths via GUI
selpath = uigetdir;
a = folder_name{1,1};

if ~exist(a, 'dir')
    outputpath = fullfile(selpath, folder_name{1});
    [status, msg] = mkdir(outputpath);
    cd(fullfile(selpath, folder_name{1}));
elseif exist(a, 'dir')
    while exist(a, 'dir')
        folder_name_warning = strcat('A folder with that name already e',...
            'xists,please choose another name');
        folder_name = inputdlg(folder_name_warning);
        a = folder_name{1,1};
    end
    outputpath = fullfile(selpath, folder_name{1});
    [status, msg] = mkdir(outputpath);
    cd(fullfile(selpath, folder_name{1}))
end


%% #An optional 3D viewer for the images in a stack. 
% Activate Code to activate 3D viewer
%imtool3D must be opened and run in order to function properly
%https://ch.mathworks.com/matlabcentral/fileexchange/40753-imtool3d
% 
% title_1 = 'Coordinates';
% message = 'please use the following 3D viewer to determine 
%the Coordinates of the spherical ROI. The coordinates of the mouse 
%location are displayed in the bottom left corner';
% f = msgbox(message, title_1);
% 
% myFolder = selpath;
% filePattern = fullfile(myFolder, '*.png');
% pngFiles = dir(filePattern);
% for k = 1:length(pngFiles)
%   baseFileName = pngFiles(k).name;
%   fullFileName = fullfile(myFolder, baseFileName);
%   fprintf(1, 'Now reading %s\n', fullFileName);
%   image = imread(fullFileName);
%   imshow(image);  % Display image.
%   imageArray(:,:,k) = image;
%   drawnow; % Force display to update immediately.
% end
% close;
% tool = imtool3D(imageArray);
% position = getPosition(tool);

%% Entering the Radius of the VOI via GUI

answer = questdlg('Do you want to define the VOI by voxels or mm?', ...
	'Radius', ...
	'Voxel','mm', 'Cancel');
% Handle response
switch answer
    case 'Voxel'
       prompt_4 = {'Enter sphere radius in voxel'};
       title_4 = 'VOI Radius';
       radius_sphere = inputdlg(prompt_4,title_4);
       radius_sphere = str2double(radius_sphere{1,1});
    case 'mm'
        prompt_2 = 'Enter scan resolution in µm';
        title_2 = 'Scan Resolution';
        dims = [1 35];
        def_input = {'2.79'};
        scan_res = inputdlg(prompt_2,title_2,dims,def_input);
        prompt_4 = {'Enter sphere radius in mm'};
        title_4 = 'VOI Z';
        radius_sphere = inputdlg(prompt_4,title_4);
        radius = ((str2double(radius_sphere{1,1}))*1000)/ ...
            str2double(scan_res{1,1});
        radius_sphere = round(radius, 0);
        clear radius;
end

message = horzcat('Enter your radius manually into the script by', ...
    'defining radius_sphere as variable (double) and manually advance the',...
    'script');
if isempty(answer) == 1
    msgbox(message);
    return
end

%% Defining the Z coordinate via GUI if not already known
% Options are by lower bound of sphere, by upper bound thereof or relative
% to two points of interest, meaning z.B. 0.7 the distance from lower bound
% to upper bound of e.g. a heart chamber or similar

answer = questdlg('Do you have the Z coordinate of the centre already?',...
    'Yes', 'No');
% Handle response
switch answer
    case 'No'
        Z_known = 0;
      question = horzcat('Do you want to define the VOI by the lower bound,',...
            'upper bound or relative to two points?');
            answer_1 = questdlg(...
                question, 'Z Coordinate', ...
                'Upper bound','Lower bound','2 Points of Interest', 'cancel');
            switch answer_1
                case 'Upper bound'
                    prompt_1 = {'Enter upper bound of sphere'};
                    title_2 = 'VOI Z';
                    Z_upper = inputdlg(prompt_1,title_2);
                    Z_centre = ((str2double(Z_upper{1,1}))- ...
                        (str2double(radius_sphere{1,1})));
                    clear Z_upper;
                case 'Lower bound'
                    prompt_1 = {'Enter lower bound of sphere'};
                    title_2 = 'VOI Z';
                    Z_lower = inputdlg(prompt_1,title_2);
                    Z_centre = ((str2double(Z_lower{1,1}))+(...
                        str2double(radius_sphere{1,1})));
                    clear Z_lower;
                case '2 Points of Interest'
                    prompt_6 = {'Enter Z of lower Point',...
                        'Enter Z of upper point', ...
                        'Enter relative distance between the points as value between 0 and 1'};
                    title_6 = 'VOI Z';
                    Z_rel = inputdlg(prompt_6,title_6);
                    Z_centre = (((str2double(Z_rel{2,1}))-(str2double(Z_rel{1,1})))* ...
                        (str2double(Z_rel{3,1})))+(str2double(Z_rel{1,1}));
            end
            message = horzcat('Enter your coordinate manually into the script by', ...
                'defining Z_centre as variable (double) and manually advance the',...
                'script');
            if isempty(answer_1) == 1
                msgbox(message);
                return
            end
    case 'Yes'
        Z_known = 1;
end

message = horzcat('Enter your radius manually into the script by', ...
    'defining radius_sphere as variable (double) and manually advance the',...
    'script');
if isempty(answer) == 1
    msgbox(message);
    return
end

%% Entering the X / Y Coordinates of the Center of the VOI via GUI

% Fiji coordinates --> Data viewer / normal coordinates
coordinate_means_fiji = 0;
if exist('counter')==0
    question = horzcat('Did you obtain your Coordinates with Fiji or Data Viewer?');
    answer_7 = questdlg( ...
        question, ...
        'Coordinate', ...
        'Fiji','Data Viewer', 'cancel');
    % Handle response
    switch answer_7
        case 'Fiji'
            coordinate_means_fiji = 1;
            prompt_res = {'Enter width of image', 'Enter height of image'};
            title_res = 'width_height';
            Resolution = inputdlg(prompt_res,title_res);
        case 'Data Viewer'
            coordinate_means_fiji = 0;
    end

    message = horzcat('Enter your coordinate manually into the script by', ...
        'defining Z_centre as variable (double) and manually advance the',...
        'script');
    if isempty(answer_7) == 1
        msgbox(message);
        return
    end
end


if coordinate_means_fiji == 1
    if Z_known == 1
        prompt_coord_fiji = {'Enter X', 'Enter Y','Enter Z'};
        title_coord_fiji = 'Fiji Coordinates';
        Coordinates_Fiji = inputdlg(prompt_coord_fiji,title_coord_fiji);
        Coordinates{1} = Coordinates_Fiji{1};
        Coordinates{2} = Resolution{2} - Coordinates_Fiji{2};
        clear Coordinates_Fiji;
        Z_centre = (str2num(Coordinates{3, 1}));
    else 
        prompt_3 = {'Enter X', 'Enter Y'};
        title_3 = 'VOI Coordinates';
        Coordinates = inputdlg(prompt_3,title_3);
        Coordinates{1} = Coordinates_Fiji{1};
        Coordinates{2} = Resolution{2} - Coordinates_Fiji{2};
        clear Coordinates_Fiji;
    end
else
    if Z_known == 1
        prompt_3 = {'Enter X', 'Enter Y', 'Enter Z'};
        title_3 = 'VOI Coordinates';
        Coordinates = inputdlg(prompt_3,title_3);
        Z_centre = (str2num(Coordinates{3, 1}));
    else 
        prompt_3 = {'Enter X', 'Enter Y'};
        title_3 = 'VOI Coordinates';
        Coordinates = inputdlg(prompt_3,title_3);    
    end
end




%% Get list of all BMP files in the specified directory

% DIR returns as a structure array.  You will need to use () and . to get
% the file names.

imagefiles = dir(fullfile(selpath, '*.png'));      
nfiles = length(imagefiles);    % Number of files found

%% Take the last 6 digits of the filename of the first file in the 
% folder without the file name extension
% define those as lowest coordinate and move the coordinate system
% accordingly

name = imagefiles(1).name;
name_length = length(name);
Z_min = ('000000');

for i = (name_length-10):(name_length-4)
    if isnan(str2num(name(i))) == 0
        Z_min(i-name_length+11) = (name(i));
    end
end

Z_min = str2double(Z_min);
Z_centre = Z_centre-Z_min;
clear name, clear name_length;

Z_max = Z_min+nfiles-1;

%% Failsafes so that the user is notified if the sphere would exceed Stack 
% dimensions
message_X_high = strcat('the defined sphere will be too far in X. ',...
    'Please revise X coordinate and / or radius');
message_X_low = strcat('the defined sphere will be too low in X. ',...
    'Please revise X coordinate and / or radius');
message_Y_high = strcat('the defined sphere will be too far in Y. ',...
    'Please revise X coordinate and / or radius');
message_Y_low = strcat('the defined sphere will be too low in Y. ',...
    'Please revise X coordinate and / or radius');
message_z_high = strcat('the defined sphere will be too high in Z. ',...
    'Please revise Z coordinate and / or radius');
message_z_low = strcat('the defined sphere will be too low in Z. ',...
    'Please revise Z coordinate and / or radius');


%X failsafes
a=imagefiles(12).name;
some_image = imread(fullfile(selpath,a), 'png');
size_some_image = size(some_image);

if  (str2num(Coordinates{1,1}) + radius_sphere) > size_some_image(1)
    uiwait(msgbox(message_X_high));
        prompt_X_high = {'Enter Centre in X', 'Enter radius'};
        title_X_high = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Coordinates{1,1});
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_X_and_radius = inputdlg(prompt_X_high,title_X_high,dims,def_input);
        Coordinates{1,1} = (str2num(new_X_and_radius{1, 1}));
        radius_sphere = (str2num(new_X_and_radius{2, 1}));
end

if (str2num(Coordinates{1,1}) - radius_sphere) < 0
    uiwait(msgboX(message_X_low));
        prompt_X_low = {'Enter Centre in X', 'Enter radius'};
        title_X_low = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Coordinates{1,1});
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_X_and_radius = inputdlg(prompt_X_low,title_X_low,dims,def_input);
        Coordinates{1,1} = (str2num(new_X_and_radius{1, 1}));
        radius_sphere = (str2num(new_X_and_radius{2, 1}));
end

% Y failsafes

if  (str2num(Coordinates{2,1}) + radius_sphere) > size_some_image(2)
    uiwait(msgbox(message_Y_high));
        prompt_Y_high = {'Enter Centre in Y', 'Enter radius'};
        title_Y_high = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Coordinates{1,1});
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_Y_and_radius = inputdlg(prompt_Y_high,title_Y_high,dims,def_input);
        Coordinates{2,1} = (str2num(new_Y_and_radius{1, 1}));
        radius_sphere = (str2num(new_Y_and_radius{2, 1}));
end

if  (str2num(Coordinates{2,1}) - radius_sphere) < 0
    uiwait(msgbox(message_Y_low));
        prompt_Y_low = {'Enter Centre in Y', 'Enter radius'};
        title_Y_low = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Coordinates{1,1});
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_Y_and_radius = inputdlg(prompt_Y_low,title_Y_low,dims,def_input);
        Coordinates{2,1} = (str2num(new_Y_and_radius{1, 1}));
        radius_sphere = (str2num(new_Y_and_radius{2, 1}));
end

%Z failsafes

if Z_max <(Z_centre+radius_sphere)
    uiwait(msgbox(message_z_high));
        prompt_Z_high = {'Enter Centre in Z', 'Enter radius'};
        title_Z_high = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Z_centre);
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_Z_and_radius = inputdlg(prompt_Z_high,title_Z_high,dims,def_input);
        Z_centre = (str2num(new_Z_and_radius{1, 1}));
        radius_sphere = (str2num(new_Z_and_radius{2, 1}));
end

if Z_min >(Z_centre-radius_sphere)
    uiwait(msgbox(message_z_low));
        prompt_Z_low = {'Enter Centre in Z', 'Enter radius'};
        title_Z_low = 'VOI Coordinates';
        dims = [1 35];
        a = num2str(Z_centre);
        b = num2str(radius_sphere);
        def_input = {a,b};
        new_Z_and_radius = inputdlg(prompt_Z_low,title_Z_low,dims,def_input);
        Z_centre = (str2num(new_Z_and_radius{1, 1}));
        radius_sphere = (str2num(new_Z_and_radius{2, 1}));
end

clear dims, clear a, clear b;

%% start stopwatch
tic

%% loading the files relevant to the VOI

f = waitbar(0,'Loading Images, please wait'); %a progress bar
%load only the relevant files
for ii=(Z_centre-radius_sphere):(Z_centre+radius_sphere) 
   waitbar(((ii-Z_centre+radius_sphere)/(2*radius_sphere)));
   currentfilename = imagefiles(ii).name;
   currentimage = imread(fullfile(selpath,currentfilename), 'png');
   images{ii} = currentimage;
end
images(cellfun('isempty', images)) = [];
close(f)

%% create a vector containing the radii for cropping

%evenly spaced vector for input in semi circle function
x = linspace(0,(2*radius_sphere),(2*radius_sphere)+1);
%semi circle function vor the respective radii
radii = ((radius_sphere^2-(x-radius_sphere).^2).^0.5);  

%% construct a cell array with all the relevant information 
% for image processing

radii_cell = num2cell(radii);
imagefile_ca = struct2cell(imagefiles);
imagefile_ca = imagefile_ca(1,(Z_centre-radius_sphere):...
    (Z_centre+radius_sphere));
ca_crop = [images; radii_cell; imagefile_ca(1,:)];
clear radii_cell, clear images, clear radii, clear x, clear imagfile_ca;

%% define cropping variables 

imageSize = size(currentimage);
X_min = ((str2num(Coordinates{1, 1}))-(radius_sphere));
Y_min = (length(currentimage)-(str2num(Coordinates{2, 1}))-(radius_sphere));
width = 2*radius_sphere;
height = 2*radius_sphere;

%% get the pictures, crop them and save them in the specified folder
f = waitbar(0,'Cropping images, please wait'); %a progress bar

for i = 1:length(ca_crop(1,:))
    waitbar((i/length(ca_crop)));
    I = imread(fullfile(selpath,ca_crop{3,i}), 'png'); %read images
    % crop images squared, so that only necessary part remains 
    % (this enables faster processing)
    I2 = imcrop(I,[X_min, Y_min, width, height]); 
    imageSize_2 = size(I2);
    % center and radius of circle ([c_row, c_col, r])
    ci = [imageSize_2(1)/2, imageSize_2(1)/2, ca_crop{2,i}];     
    % define Coordinate system based on the Coordinates specified
    [xx,yy] = ndgrid((1:imageSize_2(1))-ci(1),(1:imageSize_2(2))-ci(2));
    % create mask for cropping. Unused pixels will be zero
    mask = uint8((xx.^2 + yy.^2)<ci(3)^2); 
    croppedImage = uint8(zeros(size(I2)));
    croppedImage(:,:,1) = I2(:,:,1).*mask; %applying mask to image
    imshow(croppedImage); % show the cropped image 
    % save the cropped image under the same name as the initial image ...
    % but in the newly created folder
    imwrite(croppedImage,ca_crop{3,i}, 'png', 'WriteMode','append')
    hold on
end

close(f);
close;

%% Informing the user of the successful cropping

f = msgbox((strcat(num2str((2*radius_sphere)+1), ...
    ' images have been cropped from ', {' '}, selpath, ' and saved to ', ...
    {' '},  pwd)));

%%
% return a .txt log file with the parameters for cutting the VOI

Z_centre = Z_centre+Z_min;
title_file = strcat('Logfile_', folder_name{1},'.txt');
file = fopen(title_file, 'wt');
time_elapsed = num2str(round(toc,0));
time = char(strcat({'it took'}, {' '}, {time_elapsed}, {' '},...
    {'seconds to calculate the VOI'}));

fprintf(file, '%s\n%s\n%s%s\n%s%s\n%s\n%s\n%s\n%s\n%s', title_file, ...
    '----------------------------', 'Input directory= ', ...
    selpath, 'Output directory= ', pwd, strcat('X= ', Coordinates{1, 1}),...
    strcat('Y= ', Coordinates{2, 1}), strcat('Z= ', num2str(Z_centre)), ...
    strcat('Radius= ', num2str(radius_sphere)), time);
fclose(file);

%% run times counter

if exist('counter')==1
    counter = counter +1;
else
    counter = 1;
end