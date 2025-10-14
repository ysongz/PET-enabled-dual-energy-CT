function x = proj_back_CASTOR(G, Gopt, y, numfrm)
%--------------------------------------------------------------------------
% back projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%
% Modified by Yansong Zhu to call CASToR projector for real scanner data

if nargin<4 | isempty(numfrm)
    numfrm = 1;
end

% frame-by-frame
if numfrm>1
    y = reshape(y,[length(y(:))/numfrm,numfrm]);
    
    for m = 1:numfrm
        x(:,m) = proj_back_CASTOR(G, Gopt, y(:,m), 1);
    end
    return;
end

if numfrm==1
    y = y(:);
end
if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end



if ~isfield(Gopt,'kernel')
    Gopt.kernel = [];
end

switch Gopt.mtype
    case 'matlab'
        x = G' * y;
        
    case 'EXPLORER_histogram'
        
        if(~isfield(Gopt,'geo_file'))
            geo_file='EXPLORER';
        else
            geo_file=Gopt.geo_file;
        end
        projdirection='BACKWARD';
        if(~isfield(Gopt,'imgsiz'))
            nb_vox_X=140;   
            nb_vox_Y=140;   
            nb_vox_Z=486;
        else
            nb_vox_X=Gopt.imgsiz(1);
            nb_vox_Y=Gopt.imgsiz(2);
            nb_vox_Z=Gopt.imgsiz(3);
        end
        if(~isfield(Gopt,'fovSize'))
            fovSize_X=536.067;  
            fovSize_Y=536.067;  
            fovSize_Z=1940.1;
        else
            fovSize_X=Gopt.fovSize(1);  
            fovSize_Y=Gopt.fovSize(2);  
            fovSize_Z=Gopt.fovSize(3);
        end
        if(~isfield(Gopt,'offset'))
            offset_X=0;  offset_Y=0;  offset_Z=0;
        else
           offset_X = Gopt.offset(1);   offset_Y = Gopt.offset(2);  offset_Z = Gopt.offset(3);
        end
        dataMode='Histogram'; %Histogram or ListMode
        if(Gopt.TOF_flag==0)
            TOFinputFlag='NonTOF'; %TOF or NonTOF
            nbTOFBins=1;

            % Form input histogram as Castor event format
            % one event (1 time field+ 1 data value * TOF bin + 2 crystal IDs)
            eventsize=3+nbTOFBins;
            hist_in=-1*ones(length(Gopt.timefield)*eventsize,1);
            hist_in(1:eventsize:end)=Gopt.timefield;
            hist_in(eventsize-1:eventsize:end)=Gopt.crystalID1;
            hist_in(eventsize:eventsize:end)=Gopt.crystalID2;
            pos_val=find(hist_in==-1);
            hist_in(pos_val)=y; % set init data value
            img_in=zeros(nb_vox_X,nb_vox_Y,nb_vox_Z);
            %
            x=Castor_projector_OMP(single(img_in),single(hist_in),geo_file,projdirection,...
                                nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                                fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag);
            x=double(x(:));
        else
            TOFinputFlag='TOF'; %TOF or NonTOF
            if(~isfield(Gopt,'TOFResolutionInPs'))
                TOFResolutionInPs=505;
            else
                TOFResolutionInPs=Gopt.TOFResolutionInPs;
            end
            if(~isfield(Gopt,'nbTOFBins'))
                nbTOFBins=27;
            else
                nbTOFBins=Gopt.nbTOFBins;
            end
            if(~isfield(Gopt,'TOFBinSizeInPs'))
                TOFBinSizeInPs=273;
            else
                TOFBinSizeInPs=Gopt.TOFBinSizeInPs;
            end

            % Form input histogram as Castor event format
            % one event (1 time field+ 1 data value * 23 TOF bin + 2 crystal IDs)
            eventsize=3+nbTOFBins;
            hist_in=-1*ones(length(Gopt.timefield)*eventsize,1);
            hist_in(1:eventsize:end)=Gopt.timefield;
            hist_in(eventsize-1:eventsize:end)=Gopt.crystalID1;
            hist_in(eventsize:eventsize:end)=Gopt.crystalID2;
            pos_val=find(hist_in==-1);
            hist_in(pos_val)=y; % set init data value
            img_in=zeros(nb_vox_X,nb_vox_Y,nb_vox_Z);
            %

            x=Castor_projector_OMP(single(img_in),single(hist_in),geo_file,projdirection,...
                                nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                                fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag,TOFResolutionInPs,...
                                nbTOFBins,TOFBinSizeInPs);
            x=double(x(:));
        end 
        
end   

if not(isempty(Gopt.kernel))
    x = Gopt.kernel' * x;
end
x=single(x);
