% createiristemplate - generates a biometric template from an iris in
% an eye image.
%
% Usage: 
% [template, mask] = createiristemplate(eyeimage_filename)
%
% Arguments:
%	eyeimage_filename   - the file name of the eye image
%
% Output:
%	template		    - the binary iris biometric template
%	mask			    - the binary iris noise mask
%
% Author: 
% Libor Masek
% masekl01@csse.uwa.edu.au
% School of Computer Science & Software Engineering
% The University of Western Australia
% November 2003

function [template, mask, conf] = createiristemplate(eyeimage_filename, repo)

% path for writing diagnostic images
global DIAGPATH
DIAGPATH = 'diagnostics\';

%normalisation parameters
radial_res = 20;
angular_res = 240;
% with these settings a 9600 bit iris template is
% created

%feature encoding parameters
nscales=1;
minWaveLength=18;
mult=1; % not applicable if using nscales = 1
sigmaOnf=0.5;

eyeimage = imread([repo eyeimage_filename]); 

savefile = ['mat_data/', eyeimage_filename,'-houghpara.mat'];
[stat,mess]=fileattrib(savefile);

if stat == 1
    % if this file has been processed before
    % then load the circle parameters and
    % noise information for that file.
    load(savefile);
    %disp(['load ' savefile]);
else
    
    % if this file has not been processed before
    % then perform automatic segmentation and
    % save the results to a file
    
    [circleiris circlepupil imagewithnoise] = segmentiris(eyeimage);
    save(savefile,'circleiris','circlepupil','imagewithnoise');
    return; 
end

% WRITE NOISE IMAGE
%

imagewithnoise2 = uint8(imagewithnoise);
imagewithcircles = uint8(eyeimage);

%get pixel coords for circle around iris
[x,y] = circlecoords([circleiris(2),circleiris(1)],circleiris(3),size(eyeimage));
ind2 = sub2ind(size(eyeimage),double(y),double(x)); 

%get pixel coords for circle around pupil
[xp,yp] = circlecoords([circlepupil(2),circlepupil(1)],circlepupil(3),size(eyeimage));
ind1 = sub2ind(size(eyeimage),double(yp),double(xp));


% Write noise regions
imagewithnoise2(ind2) = 255;
imagewithnoise2(ind1) = 255;
% Write circles overlayed
imagewithcircles(ind2) = 255;
imagewithcircles(ind1) = 255;
imwrite(imagewithnoise2,[DIAGPATH, eyeimage_filename, '-noise.jpg'],'jpg');
imwrite(imagewithcircles,[DIAGPATH, eyeimage_filename, '-segmented.jpg'],'jpg');

% perform normalisation

[polar_array noise_array] = normaliseiris(imagewithnoise, circleiris(2),...
    circleiris(1), circleiris(3), circlepupil(2), circlepupil(1), circlepupil(3),eyeimage_filename, radial_res, angular_res);

% WRITE NORMALISED PATTERN, AND NOISE PATTERN
imwrite(polar_array,[DIAGPATH, eyeimage_filename, '-polar.jpg'],'jpg');
imwrite(noise_array,[DIAGPATH, eyeimage_filename, '-polarnoise.jpg'],'jpg');

% perform feature encoding
[template, mask, conf] = encode(polar_array, noise_array, nscales, minWaveLength, mult, sigmaOnf); 