function [MSEmin,KaOPT,bOPT]=multi_global_kb2(data_dir,extens,refsp,rm,rb,sigLB,sigDK,numBS,Rav,bMIN,bMAX,KaMIN,KaMAX,iter)

% original version: 4/7/2010
% revised: 10/25/2021
% R. Dyche Mullins
% please cite: Hansen SD, Mullins RD (2010) VASP is a processive actin 
% polymerase that requires monomeric actin for barbed end association. 
% J. Cell Biol. 191(3):571–584.

% This function performs a global grid search to find a least-squares, best
% fit solution for heterologous, multi-protein interactions in the 
% analytical ultracentrifuge. The function assumes that a labeled molecule 
% binds to multiple, unlabled (hence invisible) molecules. Additionally the
% function assumes that data have been collected at three different
% centrifuge speeds in a six-channel centerpiece with three sample wells. 
% The concentration of the dark molecule can vary but the unlabeled
% molecule should be present at the same concentration in all channels. 
% Finally, to reduce the number of search parameters, the model assumes
% that the cooperative (or anti-cooperative) effect is the same for every
% molecule that binds after the first one (i.e. 2nd, 3rd, 4th, 5th, etc.).

% OUTPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSEmin -- minimum value of mean square error measured over search grid
% KaOPT -- association equilibrium constant at minimum mean square error
% bOPT -- cooperativity constant at minumum mean square error

% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_dir -- name of the directory containing the data sets
% extens -- text string containing the file extension to open (e.g. 'RA1')
% refsp -- reference speed (RPM) for all the data analysis
% rm -- vector containing meniscus positions of all data sets
% rb -- vector containing base positions of all data sets
% sigLB -- sigma factor for the labeled species
% sigdk -- sigma factor for the dark (unlabeled) species
% numBS -- number of binding sites for the dark species
% Rav -- molar ratio of dark to labeled protein in each dataset
% bMIN and bMAX -- min and max values of cooperativity search axis
% KaMIN and KaMAX -- min and max values of Ka search axis
% iter -- number of points along each axis in the 2D grid search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sig -- a vector containing all the sigma factors for fitting
sig = zeros(numBS+1,1);
sig(1) = sigLB;
for i=2:numBS+1
    % add the bound dark molecules sequentially
    sig(i)=sig(i-1)+sigDK;
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sig=[1.779; 2.351; 2.923; 3.495; 4.068];
%Rav=[1.4; 1.4; 1.4; 17.5; 17.5; 17.5];   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of datasets and file names from the data directory
[dfiles,fcontent,nds] = get_datafile_names(data_dir,extens);
% nds is the number of datasets to analyze
% dfiles is a cell array holding the names of the files to analyze
% write out the data file names to the command window
dfiles  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sig is a vector containing the sigma factors for fitting
sig = zeros(numBS+1,1);
sig(1) = sigLB;
for i=2:numBS+1
    % add the bound dark molecules sequentially
    sig(i)=sig(i-1)+sigDK;
end   
% write out the sigma values to the command window
sig

% create a cell array to hold all the data and parameters
dset=cell(nds,9);
% dset{*,1} is an mx2 matrix containing the x and y values of the dataset
% dset{*,2} is the number of points in the dataset (scalar)
% dset{*,3} is the position of the meniscus (scalar)
% dset{*,4} is the position of the base of the cell (scalar)
% dset{*,5} is the speed factor associated with the data
% dset{*,6} is the area under the curve (rm to rb - before normalization)
% dset{*,7} are expected values from the linear fitting
% dset{*,8} ae the residuals associated with the fit
% dset{*,9} are the ratios of actin to VASP for each data set

% create a cell array to hold the raw data and final fit
rawdat=cell(nds,1);
% rawdat{1} is an mx4 matrix containing raw x, y, fitted y, and residuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   First we will load all of the data with a given file extension
%   found in a given directory and display the raw data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the datasets
for i=1:nds
    % open a dataset and read out the data and headers
    [x,y,err,speed,temp,sh1,sh2] = getdata(dfiles{i},data_dir);
    rawdat{i}(:,1)=x;   % do we still need rawdat?
    rawdat{i}(:,2)=y;   % do we still need rawdat?
    dset{i,1}(:,1)=x;   
    dset{i,1}(:,2)=y;
    npts = length(x);
    dset{i,2}=npts;
%    dset{i,3}=rawdat{i}(1,1);     % meniscus set to first point
%    dset{i,4}=rawdat{i}(npts,1);  % base set to last point
    dset{i,3}=rm(i);              % position of the meniscus from cmd line
    dset{i,4}=rb(i);              % position of the base from cmd line
    dset{i,5}=speed*speed/(refsp*refsp);    % speed factor for data set
    dset{i,9}=Rav(i);             % for now these come from command line
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrapolate the data curves out to rm and rb and then calculate the areas
% under the expanded data curves. After that, normalize the data so that
% each curve covers a total area of 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: careful about using 'npts' as the number of points in a dataset.
% Remember that it is a local placeholder and should only be used inside 
% the loops.
dum=cell(nds,1);    % create a cell to hold expanded version of the data
for i=1:nds
    npts=dset{i,2}+2;   % add two points to the data set: one at the
                        % meniscus and another at the base of the channel
    dum{i}(:,:)=zeros(npts,2);
    dum{i}(1,1)=dset{i,3};      % set x-value of meniscus position
    dum{i}(npts,1)=dset{i,4};   % set x-value of base position
    
    % fill most of the dummy array with data
    for j=1:npts-2
        dum{i}(j+1,1)=dset{i,1}(j,1);
        dum{i}(j+1,2)=dset{i,1}(j,2);
    end

    % nxpl is the number of points to use in the quadratic extrapolation
    % xwin and ywin are vectors to hold the points for fitting
    nxpl=16;
    xwin=zeros(nxpl,1);
    ywin=zeros(nxpl,1);
    % use linear extrapolation to find y value of the meniscus
    for j=1:nxpl
        xwin(j)=dset{i,1}(j,1);
        ywin(j)=dset{i,1}(j,2);
    end
    cdum=polyfit(xwin,ywin,1);        % do linear fit of first nxpl points
    dum{i}(1,2)=cdum(1)*dum{i}(1,1)+cdum(2); % calculate meniscus from fit
    if(dum{i}(1,2)<0)
        dum{i}(1,2)=0;
    end
        
    % use quadratic extrapolate to the y value at the base
    for j=1:nxpl
        xwin(j)=dset{i,1}(dset{i,2}-nxpl+j,1);
        ywin(j)=dset{i,1}(dset{i,2}-nxpl+j,2);
    end
    % do linear fit of first nxpl points
    cdum=polyfit(xwin,ywin,2);
    % use the fit to find the base of the cell
    dum{i}(npts,2)=cdum(1)*dum{i}(npts,1)*dum{i}(npts,1)+cdum(2)*dum{i}(npts,1)+cdum(3);
    
    % after all this effort to find points on the boundary, calculate the
    % areas under the expanded data curves.
    area=0;
    area=trapz(dum{i}(:,1),dum{i}(:,2));
    area
    dset{i,6}=area;
    dset{i,1}(:,2)=dset{i,1}(:,2)/area; % normalize the y data by the area
    dum{i}(:,2)=dum{i}(:,2)/area;       % we will use this later
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at the data to make sure that each dataset was given the correct
% meniscus and base positions, and that extrapolation didn't do anything
% pathological
figure(1);
hold on;
for i=1:nds
    plot(dum{i}(:,1),dum{i}(:,2),'blue');
    plot(dset{i,1}(:,1),dset{i,1}(:,2),'red');
end
sttl = 'data sets with extrapolated boundaries';
    title(sttl);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of data expansion, area calculation, and normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize area constants for fitting calculations. I(j,k) is the kth
% species in the ith dataset fit. REMEMBER: this calcualtion is performed 
% in real space, radial coordinates (not r^2/2 coordinates). It calls the 
% function int_radial_exp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iv=zeros(numBS+1,nds);
Ia=zeros(numBS+1,nds);
% calculate the areas under the curve for each molecular species (k) in 
% each data set (i). 
for i=1:nds
   for k=1:numBS+1      
      Iv(k,i)=int_radial_exp(1,dset{i,5}*sig(k),dset{i,3},dset{i,4});
      if(k==1)
          % the first component of the dark species is sigDK
          % currently this is the only Ia value used in the code
          Ia(k,i)=int_radial_exp(1,dset{i,5}*sigDK,dset{i,3},dset{i,4});
      else
          Ia(k,i)=(k-1)*int_radial_exp(1,dset{i,5}*sig(k),dset{i,3},dset{i,4});
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform the x-axis to linearize the data. This makes some of the 
% mathematics prettier but is not essential for the algorithm to work. 
totpts=0;
for i=1:nds
    % calculate the total number of points in all datasets
    totpts=totpts+dset{i,2};
    for j=1:dset{i,2}
        dset{i,1}(j,1)=dset{i,1}(j,1)*dset{i,1}(j,1)/2 - dset{i,3}*dset{i,3}/2;
    end
end

% initialize vectors to hold the apmlitudes of the fitting functions
ao=zeros(nds,1);
vo=zeros(nds,1);

% calculate the deltas for the Ka (dKa) x beta (db) search
dKa=(KaMAX-KaMIN)/iter;
db=(bMAX-bMIN)/iter;

% MSE is a vector holding the sum of the squared residuals 
MSE=zeros(iter+1,iter+1);
MSEgrid=zeros(iter+1,2);
% load up the x,y values for the 3D MSE plot
for l=1:iter+1
    MSEgrid(l,1)=KaMIN+(l-1)*dKa;
    MSEgrid(l,2)=bMIN+(l-1)*db;
end

% kn -- exponents for the cooperativitiy factor, used in the fitting fxn
kn=[0; 0; 1; 3; 6; 10; 15; 21; 28];
% Now do a 2-D grid search over the fitting parameters
for m=1:iter+1
    % set a new cooperativity constant
    beta=bMIN+(m-1)*db;
    for l=1:iter+1
    
        % set a new affinity constant
        Ka1=KaMIN+(l-1)*dKa;
    
        % now build the constants for calculating Vo and Ao as a function 
        % of Ka1 and then do the calculation. This is turns out to be the 
        % solution of a numBS+1 order polynomial equation, which we solve 
        % numerically. Next, compute the expected values and the residuals 
        % of the fit
        C=zeros(nds,numBS+2);
        R=zeros(nds,numBS+1);
        for i=1:nds
            % calculate the coefficients of the numBS+1 order polynomial 
            % equation that describes the amplitudes of actin-containing 
            % components
            C(i,numBS+2)=-Rav(i)*Iv(1,i);
            C(i,1)=(beta^kn(numBS+1))*(Ka1^numBS)*Ia(1,i)*Iv(numBS+1,i);
            for k=1:numBS
                C(i,numBS+2-k)=(k-Rav(i))*(beta^kn(k+1))*(Ka1^k)*Iv(k+1,i) + (beta^kn(k))*(Ka1^(k-1))*Ia(1,i)*Iv(k,i);
            end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%% calculate the amplitude coefficients - ao(i) %%%%%%%%%%%%%%%%%
            R(i,:)=roots(C(i,:));
            % look for the real positive roots
            for k=1:numBS+1
                if(imag(R(i,k))==0 && real(R(i,k))>0)
                    ao(i)=R(i,k);
                end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % calculate vo from the ao value
            voIMsum=Iv(1,i);
            for qi=1:numBS
                voIMsum = voIMsum + Iv(qi+1,i)*(beta^kn(qi+1))*(ao(i)*Ka1)^qi;
            end
            vo(i)=1/voIMsum;
    
            % use the newly calculated coefficients to compute the residuals
            dset{i,7}=fit_fxn(dset{i,1}(:,1),ao(i),vo(i),dset{i,5}*sig,Ka1,beta,numBS,kn);
            dset{i,8}=dset{i,1}(:,2)-dset{i,7};
    
            % compute the new mean square error (MSE)
            MSE(l,m)=MSE(l,m)+dset{i,8}'*dset{i,8};
        end
        if(m==1&&l==1)
            MSEmin=MSE(1,1);
            mOPT=m;
            lOPT=l;
        else
            if(MSE(l,m)<MSEmin)
                MSEmin=MSE(l,m);
                mOPT=m;
                lOPT=l;
            end
        end
        MSE(l,m)=log(MSE(l,m));
    end
end

% Now go back and generate fitted data and residuals for the optimal fit
% parameters
bOPT=MSEgrid(mOPT,2);   
KaOPT=MSEgrid(lOPT,1);
% Vtest and Atest are the areas under the fitted curves, which should be
% 1.0 for all of the V fits and Rav(i) for all of the A fits. This is a
% math check to see whether the fits make sense
Vtest=zeros(nds,1);
Atest=zeros(nds,1);
        for i=1:nds
            % calculate the coefficients of the numBS+1 order polynomial
            % equation, this time with 'optimal' parameters
            C(i,numBS+2)=-Rav(i)*Iv(1,i);
            C(i,1)=(bOPT^kn(numBS+1))*(KaOPT^numBS)*Ia(1,i)*Iv(numBS+1,i);
            kn=[0; 0; 1; 3; 6; 10; 15; 21; 28];
            for k=1:numBS
                C(i,numBS+2-k)=(k-Rav(i))*(bOPT^kn(k+1))*(KaOPT^k)*Iv(k+1,i) + (bOPT^kn(k))*(KaOPT^(k-1))*Ia(1,i)*Iv(k,i);
            end
    
            % calculate the amplitude coefficients
            R(i,:)=roots(C(i,:));
            % R(i,:)'
    
            for k=1:numBS+1
                if(imag(R(i,k))==0 && real(R(i,k))>0)
                    ao(i)=R(i,k);
                end
            end
            % calculate the Vo with optimal parameters
            voIMsum=Iv(1,i);
            for qi=1:numBS
                voIMsum = voIMsum + Iv(qi+1,i)*(bOPT^kn(qi+1))*(ao(i)*KaOPT)^qi;
            end
            vo(i)=1/voIMsum;
%             vo(i)=1/(Iv(1,i)+ao(i)*KaOPT*(Iv(2,i)+ao(i)*bOPT*KaOPT*(Iv(3,i)+ao(i)*(bOPT^2)*KaOPT*(Iv(4,i)+ao(i)*(bOPT^3)*KaOPT*Iv(5,i)))));
    
            % use the optimized coefficients to compute the residuals
            dset{i,7}=fit_fxn(dset{i,1}(:,1),ao(i),vo(i),dset{i,5}*sig,KaOPT,bOPT,numBS,kn);
            dset{i,8}=dset{i,1}(:,2)-dset{i,7};  

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % check the area under the fitted curves - this is an
            % important check on the solutions of ao and vo 
%             Atest=0;
%             Vtest=0;
            for qi=1:numBS+1
                Vtest(i) = Vtest(i) + (KaOPT*ao(i))^(qi-1)*(bOPT^kn(qi))*Iv(qi,i);
                Atest(i) = Atest(i) + (qi-1)*(KaOPT*ao(i))^(qi-1)*(bOPT^kn(qi))*Iv(qi,i);
            end
            Vtest(i) = vo(i)*Vtest(i);
            Atest(i) = ao(i)*Ia(1,i) + vo(i)*Atest(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end


% write the optimal values of fitting parameters to the workspace
ao
Atest
vo
Vtest

MSEmin
KaOPT
bOPT

% now plot all the data with the fits
figure(2)
subplot(2,2,1)
hold on
for i=1:nds
    plot(rawdat{i}(:,1),dset{i,1}(:,2),'ok');
    plot(rawdat{i}(:,1),dset{i,7},'-r','LineWidth',2);
end
%axis([0.0 1.8 0.0 12.0])
sttl = 'data and fits';
    title(sttl);
xlabel('radius (cm)','FontSize',14)
ylabel('concentration (a.u.)','FontSize',14)
hold off

% plot all of the residuals below the data and fits
subplot(2,2,3);
hold on
for i=1:nds
    plot(rawdat{i}(:,1),dset{i,8},'o');
end
%axis([0.0 1.8 0.0 12.0])
sttl = 'resuduals (data-fit)';
    title(sttl);
xlabel('radius (cm)','FontSize',14)
ylabel('residuals (a.u.)','FontSize',14)
hold off

% drop the MSE surface until it touches the x-y plane
MSE=MSE-min(min(MSE));
% do a 3D plot of the MSE versus Ka and beta
subplot(2,2,2);
hold on;
contour(MSEgrid(:,2),MSEgrid(:,1),MSE,40);
surf(MSEgrid(:,2),MSEgrid(:,1),MSE);
grid on
view(190,25)
colormap hsv
sttl = 'surface plot of MSE over search grid';
    title(sttl);
xlabel('beta','FontSize',14)
ylabel('Ka','FontSize',14)
zlabel('MS Error','FontSize',14)
hold off;

subplot(2,2,4);
hold on;
contour(MSEgrid(:,2),MSEgrid(:,1),MSE,50);
grid on
colormap hsv
sttl = 'contour plot of MSE over search grid';
    title(sttl);
xlabel('beta','FontSize',14)
ylabel('aK','FontSize',14)
hold off;

return
end

function y_fxn=fit_fxn(x_fxn,Ao,Vo,sigma,Ka,bet,NBS,bexp)
    % y_fxn -- output vector
    % sigma -- reduced eff. M.W.s of the complexes present in the mixture
    % x_fxn -- input x value vector
    % Bo -- amplitude of free actin distribution
    % Vo -- amplitude of free VASP distribution
    % bet -- cooperativity coefficient
    % NBS -- number of binding sites on labeled protein
    % bexp -- vector containing exponents for cooperativity coefficient
    
    % hard-coded for now
    bexp=[0; 0; 1; 3; 6; 10; 15; 21];
    
    np=length(x_fxn);
    y_fxn=zeros(np,1);
    for m=1:np
        % the first term of the sum is the free, labeled protein
        y_fxn(m)=exp(sigma(1)*x_fxn(m));
        % now add the contributions of all the complexes
        for q=1:NBS
            y_fxn(m)=y_fxn(m) + (bet^(bexp(q+1)))*((Ao*Ka)^q)*exp(sigma(q+1)*x_fxn(m));
        end
        y_fxn(m) = Vo*y_fxn(m);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is the original version of the fit fxn, which introduced an error
    %    y_fxn(m)=Vo*(exp(sigma(1)*x_fxn(m))+Ka*Ao*(exp(sigma(2)*x_fxn(m))+bet*Ka*Ao*(exp(sigma(3)*x_fxn(m))+(bet^2)*Ka*Ao*(exp(sigma(4)*x_fxn(m))+(bet^3)*Ka*Ao*exp(sigma(5)*x_fxn(m))))));
    % below is the first corrected version of fit fxn (10/23/2021)
%         y_fxn(m)=Vo*( ...
%             exp(sigma(1)*x_fxn(m))+ ...
%             (bet^0)*((Ka*Ao)^1)*(exp(sigma(2)*x_fxn(m)))+ ...
%             (bet^1)*((Ka*Ao)^2)*(exp(sigma(3)*x_fxn(m)))+ ...
%             (bet^3)*((Ka*Ao)^3)*(exp(sigma(4)*x_fxn(m)))+ ...
%             (bet^6)*((Ka*Ao)^4)*(exp(sigma(5)*x_fxn(m))) ...
%             );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end

return
end


function IE=int_radial_exp(B,sigma,rm,rb)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % B is the amplitudes of the exponential function
    % sigma is the reduced eff. M.W.s of the components
    % rm is the radial position of the meniscus
    % rb is the radial posision of the base
    
    delta=rb-rm;
    x=rm:delta/99:rb;
    y=zeros(100,1);    % initialize the output vector
    
    for i=1:100
        y(i)=B*exp(sigma*(x(i)*x(i)/2-rm*rm/2));
    end
    
    IE=trapz(x,y);

return
end

function [namelist,fcontent,numfiles] = get_datafile_names(data_dir,extens)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get_datafile_names - goes to a directory; identifies all the .rax files
    % and loads their names into a cell string array. 
    % data_dir - where to look for the data files
    % extens - string containing the extension for the files of interest (e.g.
    % 'RA1', 'RA2', etc.
    % numfiles - is the number of files whose names are retrieved
    % fcontent - is an array containing headers and speeds
    % RDM - 10/23/2020
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % set the return directory to the current one
    ret_dir = pwd;
    % go to the specified data directory
    cd(data_dir)
    % do a linux 'ls' command and store the results in 'txtfiles'. This
    % makes txtfiles a list of strings corresponding to the file names.
    extens = ['.' extens];
    searstr = ['*' extens];    % first make a search string
    raxfiles = ls(searstr);
    
    % how many files?
    ddum = size(raxfiles);
    numfiles = ddum(1);
    
    % create an indexable list of namess that can be used to read the data
    % files or rearranged to construct new names, linked to he originals. 
    % This list of name parts gets passed on to other functions
    namelist = cell(numfiles,1);
    
    for i=1:numfiles
        namelist{i}=raxfiles(i,:);
    end
    
    % make sure raxfiles is a single row vector
    raxfiles=raxfiles';
    raxfiles=reshape(raxfiles,1,[]);
    
    % get info about the contents of the data files
    % column 1 is the file names; column 2 is the list of speeds
    fcontent = cell(numfiles,2);
    for k=1:numfiles
        % open the data file
        fileID = fopen(namelist{k},'r');
        dum = fgetl(fileID);
        fcontent{k,1} = fgetl(fileID);
        fclose(fileID);
        
        % find the speed information
        % find all the spaces in the string
        SPind = findstr(char(32),fcontent{k,1}); % ascii code for space is 32
        dumspeed = fcontent{k,1}(SPind(3)+1:SPind(4)-1);
        fcontent{k,2} = str2num(dumspeed);
        
    end
    
    cd(ret_dir)

return
end

function [x,y,err,speed,temp,sh1,sh2] = getdata(datafile,datadir)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % open an analytical ultracentrifuge data file and copy out all the data
    %
    % x - x positional data
    % y - absorbance at each position in the cell
    % err - calculated standard deviation at each point
    % speed - centrifuge speed (determined from header information
    % sh1 - first line of the header
    % sh2 - second line of the file header
    %
    % datafile - is the filename to open
    % datadir - is the directory that holds all the data files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % set the return directory to the current one
    ret_dir = pwd;
    % go to the specified data directory
    cd(datadir)
    
    % get info about the contents of the data files
    % column 1 is the file names; column 2 is the list of speeds
    
    % open the data file
    fileID = fopen(datafile,'r');
    sh1 = fgetl(fileID);    % read the first line of the header
    sh2 = fgetl(fileID);    % read the second line of the header
        
    % find the speed information
    % find all the spaces in the string
    SPind = findstr(char(32),sh2); % ascii code for space is 32
    dumparam = sh2(SPind(3)+1:SPind(4)-1);
    speed = str2num(dumparam);
    % the code above fails when the speed is below 10000 RPM so let's fix it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(speed)
        dumparam = sh2(SPind(4)+1:SPind(5)-1);
        speed = str2num(dumparam);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dumparam = sh2(SPind(2)+1:SPind(3)-1);
    temp = str2num(dumparam);
    
    % count the number of data points
    fcontent = fgetl(fileID);
    k = 0;
    while(fcontent ~= -1)
        k=k+1;
        fcontent = fgetl(fileID);
    end
    
    % set up data arrays
    numpts = k;
    x = zeros(numpts,1);
    y = zeros(numpts,1);
    err = zeros(numpts,1);
    
    fclose(fileID);
    fileID = fopen(datafile,'r');
    % skip header
    dum = fgetl(fileID);
    dum = fgetl(fileID);
    for k=1:numpts
        fcontent = fgetl(fileID);
        A = sscanf(fcontent,'  %f  %f   %f');
        x(k) = A(1);
        y(k) = A(2);
        err(k) = A(3);
    end
    
    cd(ret_dir)

return
end