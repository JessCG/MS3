% MERGEPROCESSOR    Merges Multisizer III size distributions.
%
%   Errors? Ensure you have the most recent version at https://github.com/JessCG/MS3.git
%
%   This version of the MergeProcessor incorporates lab dilutions, Gust
%   sample sizes, and the MultiSizer pump rates.  If a dilution worksheet 
%   is not detected, the MergeProcessor will operate without tracking
%   dilutions and so will not output absolute concentrations but only 
%   relative (%) concentrations.  If the dilution worksheet does exist, the
%   following distributions will be computed in terms of the suspension 
%   that existed in the Gust chamber itself:
%       Bin Diameter (um)   
%       Normalized particle Volume(%)    
%       Particle concentration(ppm)    
%       Particle Count
%
%   Input Data:
%       Multisizer data must be in .CSV files with the following naming 
%       convention:
%           CIAS2010_356783_30_2010-09-22.#m3.CSV (or .csv)
%               |      |     |       |     |
%               |      |     |       |      ---The MS3 puts this here
%               |      |     |        ---------Analysis date
%               |      |      -----------------Tube diameter (um)
%               |       -----------------------SampleID
%                ------------------------------GroupID
%
%       The data file should have 29 lines in the header.  If it does not,
%       change the NUMHEADERLINES variable in the code.
%       
%       *** Not having 29 lines could be an issue if trying to keep track
%       of dilutions since the code requires the elapsed times, as provided
%       on line 16 of the MS3 output files. To solve the issue, ensure that
%       the option "short" is NOT ticked in the data export settings of
%       your instrument. You can re-create output files by opening .#M3
%       files with the MS3 software and re-exporting data. ***
%       
%       The data files must be in the directory specified with the variable
%       DIR_DATA in the code.
%
%       The dilutions file must be in .xls format, in 'Sheet1', and should 
%       reside in the folder stipulated in the variable DIR_DILUTIONS as 
%       set in the code.
%       The dilution file must be named in the following convention:
%           EXPERIMENT_SAMPLENUMBER_YYYY-MM-DD.xls
%       For example:
%           CIAS2010_356783_2010-09-22.xls
%
%       The Multisizer flowrate file must reside in the folder specified in
%       the variable DIR_ROOT, in the code, and must be named as follows:
%           MS3_flowrates.txt
%
%   Output Data (to .csv file):
%       Bin Diam(um)   Diff Vol(%)    Vol(ppm)    Particle Count
%
%   If no dilution file is detected, the MergeProcessor will simply output 
%   the data in the first two columns, as it did in previous versions of 
%   the code.  
%
%   Required Directory Structure:
%        - Root Directory (e.g. C:\Work\Project\Merging\ ; must contain the
%             | fileMS3_flowrates.txt)
%             |  
%              ----- ToMergeProcessor\ (data files reside here)
%             |
%              ----- MergedData\ (data will be output here, in .mat and
%             |         |           .csv formats)
%             |         |
%             |          ----- figures\ (images of merged spectra are saved
%             |                             here)
%             |
%              ----- Dilutions\ (dilutions worksheets reside here)
%             |
%             |
%              ----- Matlab (recommended, but not necessary; you could keep 
%             |                 your .m files here)
%             |
%              ----- RawData (recommended, but not necessary)
%
%   Required functions:
%       The following functions must be in matlab's search path:
%       mu2phi.m, multisizer_mergebins_30_200_400.m, 
%       To edit the path, click File > Set Path.
%
%   Running the Code:
%       It is recommended that you make a copy of the latest version of the
%       merge processor (e.g. MergeProcessor_batch_v270512.m) and rename it
%       MergeProcessor.m.  Ensure that the folder in which MergeProcessor.m
%       resides is in your Matlab search path.  To edit the path, click
%       File > Set Path.
%       Change the root directory name (dir_root, in the code).
%       To run the code, type 'MergeProcessor' or 'MergeProcessor_batch' at 
%       the prompt.
%
%   OTHER NOTES:
%   This script was written to process Multisizer data from 30-, 200-, and 
%   400-micron aperture tubes.  To alter this, change the value of the 
%   variable TUBE_SIZES (and possibly the plotting axis limits).
%
%   If you wish to stray from the recommended directory structure, you will
%   need to replace the values of all 'dir_XXXX' variables in the code.
%           
%   MULTISIZER OUTPUT FORMAT:
%       data begins on line 30 (29 lines of header) 
%       Col 1. bin number
%       Col 2. lower bin diameter (microns)
%       Col 3. center bin diameter (microns)
%       Col 4. particle count
%       Col 5. normalized volume concentration (%) %
%
%
%   By John Newgard, June 2012.  





clear all
close all

warning off

dir_root = '../';

addpath([dir_root 'Matlab' filesep]);
dir_data = [dir_root 'ToMergeProcessor' filesep];
dir_merged = [dir_root 'MergedData' filesep];
dir_fig = [dir_merged 'figures' filesep];
dir_dilutions = [dir_root 'Dilutions' filesep];
fname_flowrates = 'MS3_flowrates.txt';

file_ext = '.CSV';
tube_sizes_all = [30 200 400];
numheaderlines = 29;

conv_um3_to_ml = 1e-12; % 1x10^-12 um^3 in a ml
    
% Load the Multisizer flowrate data
fnames=dir(fullfile([dir_root fname_flowrates]));
if ~isempty(fnames)
    flowdat = load([dir_root fname_flowrates]);
    I30 = find(flowdat(:,1)==30);
    I200 = find(flowdat(:,1)==200);
    I400 = find(flowdat(:,1)==400);
    MS3_flowrates = [flowdat(I30,2) flowdat(I200,2) flowdat(I400,2)];
else
    disp(['ERROR: The Multisizer flowrate data file was not found.  Ensure that the file ' fname_flowrates])
    disp(['resides in the directory ' dir_root '.'])
    return
end
%% Set figure properties
% Determine Screen Size Variable 3=Width 4=Height
SCREEN_SIZE = get(0,'ScreenSize');
SCREEN_WIDTH = SCREEN_SIZE(3);
SCREEN_HEIGHT = SCREEN_SIZE(4);
% We want the first figure, used for merging, to be 1/2 height and full width.  
% We want the others to be 1/2 of the screen width and height
FIGURE_WIDTH_BIG = fix(SCREEN_WIDTH/1.7);
FIGURE_HEIGHT_BIG = fix(SCREEN_WIDTH/1.7);
FIGURE_WIDTH = fix(SCREEN_WIDTH/2.4);
FIGURE_HEIGHT = fix(SCREEN_HEIGHT/2);

% Set the position for each figure
FIGURE_1_COORDINATES = [0.25*SCREEN_WIDTH (SCREEN_HEIGHT-FIGURE_HEIGHT-80) FIGURE_WIDTH_BIG FIGURE_HEIGHT];
FIGURE_2_COORDINATES = [0.1*SCREEN_WIDTH (SCREEN_HEIGHT-FIGURE_HEIGHT-80) FIGURE_WIDTH FIGURE_HEIGHT];
FIGURE_3_COORDINATES = [.55*SCREEN_WIDTH (SCREEN_HEIGHT-FIGURE_HEIGHT-80) FIGURE_WIDTH FIGURE_HEIGHT];

%% Find all sample names in the data folder
eval(['fnames=dir(fullfile(dir_data,''*_',num2str(tube_sizes_all(1)),'_*',file_ext,'''));']);
for ia = 1:length(fnames)
    I = find((fnames(ia).name=='_')==1);
    tempchar = fnames(ia).name;
    samplenames{ia} = tempchar(1:I(2)-1); % e.g. = 'ME24_347'
end


%% MAIN LOOP
for ia = 1:length(samplenames) % For each set of samples...
    DilutionsOK = 1;
    samplename = samplenames{ia};
    eval(['fnames=dir(fullfile(dir_data,''' samplenames{ia},'*',file_ext,'''));']); %creates a structured array containing
                    % the contents of the directory dir_data. To see what filenames are contained in FNAMES, type fnames.name
    
    fname_char = fnames.name; % e.g. ME24_347_200_2010-02-01.#M3.CSV
    temp = fname_char(length(samplename)+1:end-4);
    I = find(temp=='#');
    I2 = find(temp=='_');
    if isfinite(I)
        suffix = temp(I2(2):I-2);
    else
        suffix = temp(I2(2):end);
    end
    %eval(['exportname = ''' samplename suffix '_merged'';']); % e.g. 'ME24_347_2010-02-01_merged'
    %eval(['excelinputfname = ''' samplename suffix ';']); % e.g. 'ME24_347_2010-02-01_merged'
    exportname = [samplename suffix '_merged'];
    fname_xls = [samplename suffix '.xls'];
    
    %% Use only the number of tube sizes necessary for this sample
    if length(fnames)==2
        tube_sizes = tube_sizes_all(1:2);
        tube_str = [num2str(tube_sizes(1)) '_' num2str(tube_sizes(2))];
    else
        tube_sizes = tube_sizes_all;
        tube_str = [num2str(tube_sizes(1)) '_' num2str(tube_sizes(2)) '_' num2str(tube_sizes(3))];
    end
    ntubes = length(tube_sizes);

    %% Create input filename syntax
    fname_raw1 = fnames(2).name; % e.g. ME24_347_30_2010-02-01.CSV
    fname_raw2 = fnames(1).name; % e.g. ME24_347_200_2010-02-01.CSV
    if ntubes==3
        fname_raw3 = fnames(3).name; % e.g. ME24_347_400_2010-02-01.CSV
    end
    
    %% Check for dilutions file
    eval(['fnames=dir(fullfile(dir_dilutions,''' fname_xls '''));']);
    if isempty(fnames)
        disp('WARNING: Concentrations are being computed WITHOUT TRACKING DILUTIONS.  This is ok if you are only interested ')
        disp('           in relative size abundance in the bed and not suspended concentrations.')
        disp('         This message has appeared because a dilutions Excel file was not found in the following directory: ')
        disp(['           ' dir_dilutions])
        disp(['         The name of the missing diultions file should be: ' fname_xls])
        disp(sprintf('\n'))
        DilutionsOK = 0;
    end

    %% Read the data files and assign data to variables listed on left below
    fid = fopen([dir_data fname_raw1]);
    C = textscan(fid,'%u%f%f%u%f','Delimiter',',','HeaderLines',numheaderlines);
    fclose(fid);
    I = find(C{1}~=0);
    bnum1 = C{1}; bnum1 = bnum1(I);
    lowd1 = C{2}; lowd1 = lowd1(I);
    centd1 = C{3}; centd1 = centd1(I);
    count1 = C{4}; count1 = count1(I);
    dvol1 = C{5}; dvol1 = dvol1(I);
    
    fid = fopen([dir_data fname_raw2]);
    C = textscan(fid,'%u%f%f%u%f','Delimiter',',','HeaderLines',29);
    fclose(fid);
    I = find(C{1}~=0);
    bnum2 = C{1}; bnum2 = bnum2(I);
    lowd2 = C{2}; lowd2 = lowd2(I);
    centd2 = C{3}; centd2 = centd2(I);
    count2 = C{4}; count2 = count2(I);
    dvol2 = C{5}; dvol2 = dvol2(I);
    
    if ntubes==3
        fid = fopen([dir_data fname_raw3]);
        C = textscan(fid,'%u%f%f%u%f','Delimiter',',','HeaderLines',29);
        fclose(fid);
        I = find(C{1}~=0);
        bnum3 = C{1}; bnum3 = bnum3(I);
        lowd3 = C{2}; lowd3 = lowd3(I);
        centd3 = C{3}; centd3 = centd3(I);
        count3 = C{4}; count3 = count3(I);
        dvol3 = C{5}; dvol3 = dvol3(I);
    end
    
    %% Sort the binned data by tube diameter
    temp = [centd1(1) centd2(1)];
    if ntubes==3
        temp = [temp centd3(1)];
    end
    [junk I] = sort(temp,'ascend');
    
    if ntubes==3
        eval(['centd_mu = [centd' num2str(I(1)) ' centd' num2str(I(2)) ' centd' num2str(I(3)) '];'])
        eval(['lowd_mu = [lowd' num2str(I(1)) ' lowd' num2str(I(2)) ' lowd' num2str(I(3)) '];'])
        eval(['dvol = [dvol' num2str(I(1)) ' dvol' num2str(I(2)) ' dvol' num2str(I(3)) '];'])
        eval(['numparticles = [count' num2str(I(1)) ' count' num2str(I(2)) ' count' num2str(I(3)) '];'])
    else
        eval(['centd_mu = [centd' num2str(I(1)) ' centd' num2str(I(2)) '];'])
        eval(['lowd_mu = [lowd' num2str(I(1)) ' lowd' num2str(I(2)) '];'])
        eval(['dvol = [dvol' num2str(I(1)) ' dvol' num2str(I(2)) '];'])
        eval(['numparticles = [count' num2str(I(1)) ' count' num2str(I(2)) '];'])
    end
    centd_phi = mu2phi(centd_mu);
    lowd_phi = mu2phi(lowd_mu);
    
    % Get merge bins
    centd_mu_new = multisizer_mergebins_30_200_400;
    %eval(['[centd_mu_new,lowerd_mu,upperd_mu] = multisizer_mergebins_' tube_str ';'])
    %eval(['centd_mu_new = multisizer_mergebins_' tube_str ';'])
    
    
    % Get the volumes and counts that correspond to the merge bins.
    %   - merged volumes correpond to the sum of the volumes in the 10
    %       consecutive bins centered about the merge bins
    %   - if there are less than 5 diameter bins on either side of the
    %       merge bin centre, the volume for that merged bin is given a
    %       value of NaN and the particle count is given a value of zero.
    phi = mu2phi(centd_mu_new);
    for ii = 1:ntubes
        diam = centd_phi(:,ii);
        vol = dvol(:,ii);
        num = numparticles(:,ii);
        for jj = 1:length(phi)
            [y I] = min(abs(phi(jj) - diam));
            if diam(I) < phi(jj)
                if I-4>0 & I+5<=length(diam)
                    dvol_new(jj,ii) = sum(vol(I-4:I+5));
                    numparticles_new(jj,ii) = sum(num(I-4:I+5));
                else
                    dvol_new(jj,ii) = NaN;
                    numparticles_new(jj,ii) = NaN;
                end
            else
                if I-5>0 & I+4<length(diam)
                    dvol_new(jj,ii) = sum(vol(I-5:I+4));
                    numparticles_new(jj,ii) = sum(num(I-5:I+4));
                else
                    dvol_new(jj,ii) = NaN;
                    numparticles_new(jj,ii) = NaN;
                end
            end
        end
    end
    clear diam vol num
    clear numparticles dvol
    
    % Remove spurious points with <10 counts/bin at large-diameter end of
    % each distribution.  NB. Data in diameter bins >= the smallest of
    % these bins will ALL be removed.
    for ii = 1:ntubes
        % First replace zero values with NaNs
        I = find(numparticles_new(:,ii)==0);
        numparticles_new(I,ii) = NaN;
        dvol_new(I,ii) = NaN;
        % Now find counts less than 10
        I = find(numparticles_new(:,ii)<10);
        if isfinite(I)
            I_good = [1:min(I)-1]';
        else
            I_good = [1:length(numparticles_new(:,ii))];
        end
        I_bad = find(ismember([1:length(numparticles_new)],I_good)==0);
        numparticles_new(I_bad,ii) = NaN;
        dvol_new(I_bad,ii) = NaN;
    end
    
    % Remove NaN values from distributions
%     I = find(isfinite(numparticles_new(:,1))==1);
%     dvol1 = dvol_new(I,1);
%     centd_mu_1 = centd_mu_new(I);
%     I = find(isfinite(numparticles_new(:,2))==1);
%     dvol2 = dvol_new(I,2);
%     centd_mu_2 = centd_mu_new(I);
%     if ntubes==3
%         I = find(isfinite(numparticles_new(:,3))==1);
%         dvol3 = dvol_new(I,3);
%         centd_mu_3 = centd_mu_new(I);
%         %lowerd_mu_3 = lowerd_mu(I);
%         %upperd_mu_3 = upperd_mu(I);
%     end
    
    % Replace NaN values with zeros
    dvol_new(isnan(dvol_new(:,1)),1) = 0;
    dvol1 = dvol_new(:,1);
    dvol_new(isnan(dvol_new(:,2)),2) = 0;
    dvol2 = dvol_new(:,2);
    numparticles_new(isnan(numparticles_new(:,1)),1) = 0;
    numparticles1 = numparticles_new(:,1);
    numparticles_new(isnan(numparticles_new(:,2)),2) = 0;
    numparticles2 = numparticles_new(:,2);
    if ntubes==3
        dvol_new(isnan(dvol_new(:,3)),3) = 0;
        dvol3 = dvol_new(:,3);
        numparticles_new(isnan(numparticles_new(:,3)),3) = 0;
        numparticles3 = numparticles_new(:,3);
    end
    centd1 = centd_mu_new;
    centd2 = centd_mu_new;
    if ntubes==3
        centd3 = centd_mu_new;
    end
    clear dvol_new numparticles_new
    
    
    % Get the pumped volumes for the diluted samples
    time_pumped = NaN*ones(1,3);
    for ib = 1:ntubes
        eval(['tempname = fname_raw' num2str(ib) ';'])
        fid = fopen([dir_data tempname]);
        data = textscan(fid,'%s%s','Delimiter',',');
        I = find(strcmp(data{1},'"Elapsed time:"')==1);
        if isempty(I)
            I = find(strcmp(data{1},'Elapsed time:')==1);
            temp = data{2};
            temp2 = temp{I};
            time_pumped(ib) = str2num(temp2); % [s]
        else
            temp = data{2};
            temp2 = temp{I};
            I = find(temp2=='"');
            time_pumped(ib) = str2num(temp2(I(1)+1:I(2)-1)); % [s]
        end
    end
    vol_pumped = time_pumped .* MS3_flowrates; % [ml]
    
    % Load the dilution coefficients
    if DilutionsOK
        [xls_num,xls_txt,xls_raw] = xlsread([dir_dilutions fname_xls],'Sheet1','','basic');
        Gust_vol_sampled = xls_num(2,1); %[ml]
        Gust_frac_sampled = xls_num(2,2);%[-]
        Gust_vol = Gust_vol_sampled/Gust_frac_sampled; %[ml]
        lab_vol_added = xls_num(2,3);%[ml]
        lab_dilutions = [xls_num(2,4) xls_num(2,5) xls_num(2,6)];%[-]
    end
    
    % We do not want to use %Vol Conc for merging--we want volume concentration
    % in ppm.
    % Convert nparticles to ppm assuming every particle in each size class
    % is spherical and has identically the diameter corresponding to the 
    % bin median diameter. 
    centv = 1/6*pi*centd1.^3; %um^3
    volparticles1 = numparticles1.*centv; %um^3
    volparticles2 = numparticles2.*centv;
    ppm1_dilute_ms3 = volparticles1./vol_pumped(1)*conv_um3_to_ml*1e06; 
    ppm2_dilute_ms3 = volparticles2./vol_pumped(2)*conv_um3_to_ml*1e06; 
    if ntubes==3
        volparticles3 = numparticles3.*centv;
        ppm3_dilute_ms3 = volparticles3./vol_pumped(3)*conv_um3_to_ml*1e06;
    end
    if DilutionsOK
        ppm1_stock = ppm1_dilute_ms3/lab_dilutions(1);
        ppm2_stock = ppm2_dilute_ms3/lab_dilutions(2);
        if ntubes==3            
            ppm3_stock = ppm3_dilute_ms3/lab_dilutions(3);
        end
    else
        ppm1_stock = ppm1_dilute_ms3;
        ppm2_stock = ppm2_dilute_ms3;
        if ntubes==3            
            ppm3_stock = ppm3_dilute_ms3;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    mergeok = 0;
    while mergeok == 0;
        % plot unmerged curves to determine diameters for merging
         % Set axis limits for figures
        xmin = 1; xmax = 1000;
        I1 = find(ppm1_stock>0);  I2 = find(ppm2_stock>0);  
        if ntubes==3
            I3 = find(ppm3_stock>0);
            ymin = 0.1* min([min(ppm1_stock(I1)) min(ppm2_stock(I2)) min(ppm3_stock(I3))]); ymax = 10* max(max([ppm1_stock ppm2_stock ppm3_stock]));
        else
            ymin = 0.1* min([min(ppm1_stock(I1)) min(ppm2_stock(I2))]); ymax = 10* max(max([ppm1_stock ppm2_stock]));
        end
        figure(1), clf
        set(gcf,'Units','Pixels','Position',FIGURE_1_COORDINATES);
        %loglog(centd1,dvol1,'.-r'), hold on
        %loglog(centd2,dvol2,'.-b')
        loglog(centd1,ppm1_stock,'.-r'), hold on
        loglog(centd2,ppm2_stock,'.-b')
        if ntubes==3
            %loglog(centd3,dvol3,'.-g')
            loglog(centd3,ppm3_stock,'.-g')
        end
        xlabel('Diameter (\mum)','fontsize',14)
        %ylabel('Normalized Vol. Conc. (%)','fontsize',14)
        ylabel('Vol. Conc. (ppm) in lab stock','fontsize',14)
        axis([xmin xmax ymin ymax])
        axis square
        
        %% Find bins associated with merging diameters
        if ntubes==3
            title([samplename '          Click on two suitable merge points, then hit ENTER'],'Interpreter','none','FontSize',14)
            disp('USER INPUT REQUIRED:');
            disp('  Use the crosshairs to select the x-values at which you want ');
            disp('   the red/blue, and blue/green, curves to be merged.');
            disp('  1. Line-up the vertical line (horiz line does not matter) where you want it.');
            disp('  2. Click the mouse.');
            disp('  3. Move to the 2nd merge point and click.');
            disp('  4. Hit ENTER.');
        else
            title([samplename '          Click on a suitable merge point, then hit ENTER'],'Interpreter','none','FontSize',14)
            disp('USER INPUT REQUIRED:');
            disp('  Use the crosshairs to select the x-value at which you want ');
            disp('   the red and blue curves to be merged.');
            disp('  1. Line-up the vertical line (horiz line does not matter) where you want it.');
            disp('  2. Click the mouse.');
            disp('  3. Hit ENTER.');
        end
        
        % Find the indices of the merge bins
        % Dist'ns 1 & 2
        [mergepts,junk] = ginput; % puts cross-hairs on the figure and saves the points
        if size(mergepts)~=ntubes-1
            close
            if ntubes==3
                disp('You must click on the figure TWICE.  Please re-try.')
            else
                disp('It did not work.  Try again.')
            end
            pause(3); % This pauses for 3 seconds
        else
            j1 = max(find(centd1<=mergepts(1)));
            k1 = min(find(centd2>=mergepts(1)));
            % Dist'ns 2 & 3
            if ntubes==3
                j2 = max(find(centd2<=mergepts(2)));
                k2 = min(find(centd3>=mergepts(2)));
            end
            % Close the figure
            close
            
            % MERGING SCHEME
            % 1. Scale the smallest-aperture curve to meet the middle curve and truncate it up to
            % the merge point.
            %beta1 = dvol2(k1)/dvol1(j1);
            %merge_dvol1 = beta1.*dvol1(1:j1);
            beta1 = ppm2_stock(k1)/ppm1_stock(j1);
            merge_ppm1 = beta1.*ppm1_stock(1:j1);
            merge_centd1 = centd1(1:j1);
            
            if ntubes==3
                % 2. Scale the largest-aperture curve to meet the middle curve and truncate it up to
                % the merge point.
                %beta2 = dvol2(j2)/dvol3(k2);
                %merge_dvol3 = beta2.*dvol3(k2:end);
                beta2 = ppm2_stock(j2)/ppm3_stock(k2);
                merge_ppm3 = beta2.*ppm3_stock(k2:end);
                merge_centd3 = centd3(k2:end);
            else
                beta2 = [];
                merge_ppm3 = [];
                merge_centd3 = [];
            end
            
            if ntubes==3
                % 3. Truncate the middle curve between the merge points.
                merge_ppm2 = ppm2_stock(k1:j2);
                merge_centd2 = centd2(k1:j2);
            else
                merge_ppm2 = ppm2_stock(k1:end);
                merge_centd2 = centd2(k1:end);
            end
            % 4. Merge these truncated curves
            merge_ppm = [merge_ppm1; merge_ppm2; merge_ppm3];
            merge_centd = [merge_centd1; merge_centd2; merge_centd3];
            
            % Plot the merged curve with the unmerged curves
            %         xmin = min(centd1) - 0.5*min(centd1);
            %         xmax = max(centd3) + max(centd3);
            %         ymin = min([min(merge_dvol) min(dvol3)]) - 0.5*min([min(merge_dvol) min(dvol3)]);
            %         ymax = max([dvol1; dvol2; dvol3; merge_dvol]) + max([dvol1; dvol2; dvol3; merge_dvol]);
            
            %xmin = 1; xmax = 1000; ymin = .1; ymax = 100;
            figure(2), clf
            set(gcf,'Units','Pixels','Position',FIGURE_2_COORDINATES);
            loglog(merge_centd,merge_ppm,'-k','LineWidth',2), hold on
            loglog(centd1,ppm1_stock,'-r','LineWidth',1)
            loglog(centd2,ppm2_stock,'-b','LineWidth',1)
            if ntubes==3
                loglog(centd3,ppm3_stock,'-g','LineWidth',1)
            end
            loglog([mergepts(1) mergepts(1)],[ymin ymax],'-.','Color',[.7 .7 .7])
            if ntubes==3
                loglog([mergepts(2) mergepts(2)],[ymin ymax],'-.','Color',[.7 .7 .7])
            end
            %axis([xmin xmax ymin ymax])
            axis square
            xlabel('Diameter (\mum)','fontsize',14)
            ylabel('Volume Conc. (ppm)','fontsize',14)
            if ntubes==3
                AX=legend('merged',[num2str(tube_sizes(1)) ' \mum'],[num2str(tube_sizes(2)) ' \mum']...
                    ,[num2str(tube_sizes(3)) ' \mum']);
            else
                AX=legend('merged',[num2str(tube_sizes(1)) ' \mum'],[num2str(tube_sizes(2)) ' \mum']);
            end
            set(AX,'FontSize',8), legend boxoff
            title('Accept or reject the merged curve?')
            
            
            
            %% Normalize the merged distribution
            merge_dvol = merge_ppm./sum(merge_ppm).*100; % merge ppm is relative to the stock
        
            figure(3), clf
            ymin = .01; ymax = 100;
            set(gcf,'Units','Pixels','Position',FIGURE_3_COORDINATES);
            loglog(merge_centd,merge_dvol,'-k','LineWidth',2), hold on
            loglog([mergepts(1) mergepts(1)],[ymin ymax],'-.','Color',[.7 .7 .7])
            if ntubes==3
                loglog([mergepts(2) mergepts(2)],[ymin ymax],'-.','Color',[.7 .7 .7])
            end
            axis([xmin xmax ymin ymax])
            xlabel('Diameter (\mum)','fontsize',14)
            ylabel('Normalized Vol. Conc. (%)','fontsize',14)
            title('Accept or reject the merged curve?')
            axis square
%             if min(merge_dvol) < 0.1
%                 axis([xmin xmax .01 ymax])
%                 set(gca,'PlotBoxAspectRatio',[3 4 1])
%             end
            
            answerok = 0;
            while answerok == 0;
                yesorno = input('Do you accept this merged curve (y/n)? [y]','s');
                if isempty(yesorno) | strcmp(yesorno,'y') | strcmp(yesorno,'yes') | strcmp(yesorno,'Y') | strcmp(yesorno,'Yes') | strcmp(yesorno,'YES')
                    mergeok = 1; answerok = 1;
                elseif strcmp(yesorno,'n') | strcmp(yesorno,'no') | strcmp(yesorno,'No') | strcmp(yesorno,'NO') | strcmp(yesorno,'N')
                    mergeok = 0; answerok = 1;
                else disp(['Your answer was not understood.  Please type ''y'' or ''n''.'])
                    answerok = 0;
                end
            end
            figure(2), close
            figure(3)
            % Save the figure to a jpeg file
            title(samplename,'Interpreter','none')
            print('-djpeg100','-r200',[dir_fig exportname '.jpg'])
        end
        
    end
    
    
    %% FORMAT THE MERGED DISTRIBUTIONS FOR OUTPUT
    
    % Convert merged ppm to merged numparticles
    if DilutionsOK
        merge_centv = 1/6*pi*merge_centd.^3; %um^3
        merge_vol = merge_ppm*lab_vol_added./conv_um3_to_ml/1e06; % lab_vol_added is the volume of water added to the gust sediment sample
        merge_numparticles = merge_vol./merge_centv;
    end
    merge_dvol = merge_ppm./sum(merge_ppm)*100;
    
    % Now convert stock values to Gust values
    if DilutionsOK
        gust_numparticles = merge_numparticles *Gust_vol_sampled/lab_vol_added/Gust_frac_sampled;
        gust_dvol = merge_dvol;
        gust_vol = merge_vol *Gust_vol_sampled/lab_vol_added/Gust_frac_sampled;
        gust_ppm = merge_ppm*lab_vol_added/Gust_vol_sampled;
    end
        
    %% OUTPUT THE MERGED DISTRIBUTIONS
    if DilutionsOK
        % Save text file
        data_out = [merge_centd merge_dvol gust_ppm gust_numparticles];
        dlmwrite([dir_merged exportname '.csv'],data_out,'precision',4,'delimiter',',')
        %dlmwrite([dir_merged exportname '.csv'],data_out,'precision','%15.4f','delimiter',',')
        %dlmwrite([dir_merged exportname '.txt'],data_out,'precision',4,'delimiter','\t')
        %edit([dir_merged exportname '.txt']) % This opens the output text
                                                % file in the matlab editor
        % Save .mat file
        eval(['save ' dir_merged exportname '.mat merge_centd merge_dvol gust_ppm gust_numparticles'])
    else
        data_out = [merge_centd merge_dvol];
        dlmwrite([dir_merged exportname '.csv'],data_out,'precision',4,'delimiter',',')
        %dlmwrite([dir_merged exportname '.csv'],data_out,'precision','%15.4f','delimiter',',')
        %dlmwrite([dir_merged exportname '.txt'],data_out,'precision',4,'delimiter','\t')
        %edit([dir_merged exportname '.txt']) % This opens the output text
                                                % file in the matlab editor
        % Save .mat file
        eval(['save ' dir_merged exportname '.mat merge_centd merge_dvol'])
    end
end
