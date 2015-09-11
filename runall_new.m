addpath('./', './Segmentation', './Matching', './Normal_encoding')

%only works for nscales = 1.  This is OK - this analysis appears to be
%matching (roughly) the 20% error rate observed in Libor Masek's thesis.

%PHASE 1 - run the iris data extraction - generates houghpara files for
%analysis
do_phase1 = 1;

%PHASE 2 - Collect and flatten confidence information into a single array
do_phase2 = 1;

%PHASE 3A - Fit data using Gaussians - calculate maximum number of stable
%bits using fitted values
do_phase3A = 0;

%PHASE 3B - calculate epsilon_2 by numerically analyzing data.  Uses this
%epsilon_2 to calculate maximum number of stable bits.
do_phase3B = 1;

%ENTERING PHASE 1 - EXPERIMENTS ON IRIS DATA
radial_res = 20;
angular_res = 240;

repo = '..\TemplateAging\spring_2008-lr-clean\';

if(do_phase1 == 1)
    filelist = dir([repo '*.tiff']);
    s = struct('filename','',...
        'temp',zeros(radial_res,2*angular_res),...
        'mask',zeros(radial_res,2*angular_res),...
        'conf',zeros(radial_res,2*angular_res)); 

    totallist = repmat(s,length(filelist),1);

    for i = 1:length(filelist)
        filename = filelist(i).name;
        disp(filename);
        totallist(i).filename = filename;
        try
            [totallist(i).temp, totallist(i).mask, totallist(i).conf] = createiristemplate(filename, repo); 
        end
    end

    %savefile = ['total.mat'];
    %save(savefile,'totallist');
end

%PHASE 2 - flatten data for analysis

if do_phase2 == 1
    %calculate \sigma_INTRA by calculating pixel-by-pixel differences in
    %confidence information.  Aggregate these differences into an overall list
    % This is designed to work with the CrossSensor directory structure right
    % now, where there are 6 photos each of left, right eye of person per
    % folder. This script compares image 1 to the remainder of the files in
    % the directory.
    disp('Calculating \sigma_INTRA');
    diffconfrealflat = zeros(radial_res*angular_res*(length(filelist)/2 - 1),1);
    diffconfimagflat = zeros(size(diffconfrealflat));
    diffmaskflat = zeros(size(diffconfrealflat));

    numpixelstotest = 30;
    numhistrows=5;
    pixelcol = round(rand(numpixelstotest,1) * angular_res);
    pixelrow = round(rand(numpixelstotest,1) * radial_res);
    pixelrealvals = zeros(numpixelstotest, length(filelist));
    pixelmask = totallist(1).mask(pixelrow + radial_res*pixelcol);
    pixelmask(find(pixelmask)) = NaN;
    pixelrealvals(:,1) = totallist(1).conf(pixelrow + radial_res*pixelcol) + pixelmask;

    for i = 2:length(filelist);
        %NOTE: this assumes all images are of the same eye
        [hd,shift]=gethammingdistance(...
            totallist(i).temp,totallist(i).mask,...
            totallist(1).temp, totallist(1).mask,1);

        %shift the confidence information, subtract it off
        %disp(['Hamming Distance: ' num2str(hd)]);
        sconf = shiftbits(totallist(i).conf,shift,1);
        smask = shiftbits(totallist(i).mask,shift,1);
        confdiff = totallist(1).conf - sconf;
        maskdiff = totallist(1).mask | smask;

        %store pixel values
        pixelmask = +maskdiff(pixelrow + radial_res*pixelcol);
        pixelmask(find(pixelmask)) = NaN;
        pixelrealvals(:,i) = sconf(pixelrow + radial_res*pixelcol) + pixelmask;

        %store into a flat array for analysis
        diffconfrealflat(radial_res*angular_res*(i-2)+1:...
            radial_res*angular_res*(i-1)) =...
            reshape(confdiff(:,1:2:angular_res*2),1,radial_res*angular_res);
        diffconfimagflat(radial_res*angular_res*(i-2)+1:...
            radial_res*angular_res*(i-1)) =...
            reshape(confdiff(:,2:2:angular_res*2),1,radial_res*angular_res);
        diffmaskflat(radial_res*angular_res*(i-2)+1:...
            radial_res*angular_res*(i-1)) =...
            reshape(maskdiff(:,2:2:angular_res*2),1,radial_res*angular_res);


        [hd,shift]=gethammingdistance(...
            totallist(12).temp,totallist(12).mask,...
            totallist(1).temp, totallist(1).mask,1);
        %disp(['Hamming Distance (l to r): ' num2str(hd)]);

    end


    %calculate \sigma_INTER by histogram-ing the confidence information
    disp('Calculating \sigma_INTER');
    conflistrealflat = zeros(radial_res*angular_res*length(filelist),1);
    conflistimagflat = zeros(size(conflistrealflat));
    masklistflat = zeros(size(conflistrealflat));
    masklistimagflat = zeros(size(conflistrealflat));

    for i = 0:(length(filelist)-1)
        for j = 0:radial_res-1
            conflistrealflat((angular_res*j + radial_res*angular_res*i)+1:...
                angular_res*j + angular_res*radial_res*i + angular_res) = ...    
                totallist(i+1).conf(j+1,1:2:angular_res*2);

            conflistimagflat((angular_res*j + radial_res*angular_res*i)+1:...
                angular_res*j + angular_res*radial_res*i + angular_res) = ...    
                totallist(i+1).conf(j+1,2:2:angular_res*2);

            %translate mask
            masklistflat((angular_res*j + radial_res*angular_res*i)+1:...
                angular_res*j + angular_res*radial_res*i + angular_res) = ...    
                totallist(i+1).mask(j+1,1:2:angular_res*2);
        end
    end

end

%PHASE 3A - fit data to extract sigma values.  Use this to estimate number
%of stable bits.

if do_phase3A == 1
    
    %sigma_INTER
    %Plot distributions from random pixels across 30 samples
    %figure;
    %for i = 1:numpixelstotest
    %    subplot(numhistrows, ceil(numpixelstotest/numhistrows), i);
    %    hist(pixelrealvals(i,:),-0.03:0.002:0.03);
    %end

    %plot histogram of real data
    %hist(diffconfrealflat(find(1-diffmaskrealflat)),200)
    %PESSIMISTIC
    figure;
    histfit(diffconfrealflat(find(1-diffmaskflat)),200)
    pdintra = fitdist(diffconfrealflat(find(1-diffmaskflat)),'normal');

    %fit histogram using 'fit' function
    %OPTIMISTIC
    %[y,x]=hist(diffconfrealflat(find(1-diffmaskflat)),200);
    %f1=fit(x',y','gauss1')%, 'Exclude', not(0.05>x>-0.05));
    %figure;
    %plot(f1,x,y)%,or(x<-0.05,x>0.05))

    %custom fit type for log-exp fit for gabor filtered data
    %g = fittype('K1.*exp(-log(y./K2).^2./(2*log(K3)))',...
    %   'independent',{'y'},...
    %   'coefficients',{'K1','K2','K3'});
    %f2=fit(x',y',g,'Exclude',not(0.03>x>-0.03));
    %figure;
    %plot(f2,x,y,or(x<-0.03,x>0.03))


    %3-d histgram of real and imaginary difference
    %figure;
    %hist3([diffconfrealflat(find(1-diffmaskflat)), ...
    %    diffconfimagflat(find(1-diffmaskflat))], ...
    %    'Edges',{-0.03:0.002:0.03,-0.03:0.002:0.03});
    
    %sigma_INTRA
    
    %pick out elements that are in the mask and histogram
    %figure;
    %OPTIMISTIC
    %histfit(conflistrealflat(find(1-masklistflat)),200)
    %figure;
    %histfit(conflistimagflat(find(1-masklistflat)),200)
    %pdinter = fitdist(conflistrealflat(find(1-masklistflat)),'normal');

    %fit histogram using 'fit' function
    %PESSIMISTIC
    [y,x]=hist(conflistrealflat(find(1-masklistflat)),200);
    f2=fit(x',y','gauss1')%, 'Exclude', not(0.05>x>-0.05));
    %renormalize
    f2.a1=0.0012*sum(y)/(f2.c1/sqrt(2) * sqrt(2*pi));
    figure;
    plot(f2,x,y)%,or(x<-0.05,x>0.05))
    
    
    
    %Calculate parameters based on fitted values
    sigrat = pdintra.sigma/(f2.c1/sqrt(2))
    %sigrat = 0.75;
    e2=5.3*10^(-3); %choose this value as convenient
    pstable = (1-erf(sigrat * erfinv(1 - 2*e2)))
    numbits=2*radial_res*angular_res;
    e1=5*10^(-2); %choose this value as convenient
    syms mp;
    vpasolve(numbits*pstable==mp-log(e1)+sqrt(log(e1)*(log(e1)-2*mp)))
    %mp is the number of bits available for LPN

end

%PHASE 3B - do numerical analysis to extract epsilon_2 from distribution.
%Use this epsilon_2 to compute maximum number of stable bits.

if do_phase3B == 1
    %calculate epsilon 2 from numerical data
    %give a value of p_stable a priori
    pstable = 0.04;
    %find threshold that corresponds to this stability
    [y,x]=hist(conflistrealflat(find(1-masklistflat)),400);
    npstable = 0;
    thresh = 0.1;
    while(npstable < pstable)
        thresh = thresh - (x(2)-x(1)); %decrement threshold
        %calculate the probability of stability for that threshold
        npstable = sum(y(find(abs(x) > thresh)))/sum(y);
    end
    pstable = npstable;

    %now have a threshold for p_stable (thresh)
    %calculate epsilon_2 based on probability of bit-flip
    [y2,x2]=hist(diffconfrealflat(find(1-diffmaskflat)),400);
    %multiply next line by 1/2 so we don't double-count probabilities
    e2 = sum(y2(find(abs(x2) > thresh)))/sum(y2) * 1/2
    figure;
    plot(x2,y2/(sum(y2) * (x2(2)-x2(1))));

    %Calculate number of available bits for LPN
    numbits=2*radial_res*angular_res;
    e1=5*10^(-2);
    syms mp;
    vpasolve(numbits*pstable==mp-log(e1)+sqrt(log(e1)*(log(e1)-2*mp)))
    %mp is the number of available bits for LPN
    %e2*mp is the expected number of errors (should be less than 3 for
    %the solution to be remotely practical)
   
end

%other useful code: calculating CDF from histogram:
%ycdf = zeros(length(y),1);
%for i = 1:length(y)
%    ycdf(i) = sum(y(1:i))/sum(y);
%end
%figure;
%plot(x,ycdf);