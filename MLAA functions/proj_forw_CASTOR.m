function y = proj_forw_CASTOR(G, Gopt, x, numfrm)
%--------------------------------------------------------------------------
% forward projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%
% Modified by Yansong Zhu to call CASToR projector for real scanner data

if nargin<4 | isempty(numfrm)
    numfrm = 1;
end

% frame-by-frame
if numfrm>1
    x = reshape(x,[length(x(:))/numfrm,numfrm]);
    for m = 1:numfrm
        y(:,m) = proj_forw_CASTOR(G, Gopt, x(:,m), 1);
    end
    return;
end

if numfrm==1
    x = x(:);
end
if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end

% use kernel
if ~isfield(Gopt,'kernel')
    Gopt.kernel = [];
end
if not(isempty(Gopt.kernel))
    x = Gopt.kernel * double(x);
end

switch Gopt.mtype
    case 'matlab'
        y = G * x;
        
    case 'EXPLORER_histogram'
        
        if(~isfield(Gopt,'geo_file'))
            geo_file='EXPLORER';
        else
            geo_file=Gopt.geo_file;
        end
        projdirection='FORWARD';
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
            hist_in(pos_val==-1)=0;% set init data value
            %
            x=reshape(x,nb_vox_X,nb_vox_Y,nb_vox_Z);
            hist_out=Castor_projector_OMP(single(x),single(hist_in),geo_file,projdirection,...
                               nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                               fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag);
            y=hist_out(pos_val);
            y=y(:);
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
            % one event (1 time field+ 1 data value * TOF bin + 2 crystal IDs)
            eventsize=3+nbTOFBins;
            hist_in=-1*ones(length(Gopt.timefield)*eventsize,1);
            hist_in(1:eventsize:end)=Gopt.timefield;
            hist_in(eventsize-1:eventsize:end)=Gopt.crystalID1;
            hist_in(eventsize:eventsize:end)=Gopt.crystalID2;
            pos_val=find(hist_in==-1);
            hist_in(pos_val==-1)=0; % set init data value
            %
            x=reshape(x,nb_vox_X,nb_vox_Y,nb_vox_Z);

            hist_out=Castor_projector_OMP(single(x),single(hist_in),geo_file,projdirection,...
                                nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                                fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag,TOFResolutionInPs,...
                                nbTOFBins,TOFBinSizeInPs);
            y=hist_out(pos_val);
            y=y(:);
        end
        
    case 'EXPLORER_listmode'
        
        if(~isfield(Gopt,'geo_file'))
            geo_file='EXPLORER';
        else
            geo_file=Gopt.geo_file;
        end
        projdirection='FORWARD';
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
        dataMode='ListMode'; %Histogram or ListMode
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
            hist_in(pos_val==-1)=0;% set init data value
            %
            x=reshape(x,nb_vox_X,nb_vox_Y,nb_vox_Z);
            hist_out=Castor_projector_OMP(single(x),single(hist_in),geo_file,projdirection,...
                               nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                               fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag);
            y=hist_out(pos_val);
            y=y(:);
        else
            TOFinputFlag='TOF'; %TOF or NonTOF
            if(~isfield(Gopt,'TOFResolutionInPs'))
                TOFResolutionInPs=505;
            else
                TOFResolutionInPs=Gopt.TOFResolutionInPs;
            end
            if(~isfield(Gopt,'TOFQuantizationBinSize'))
                TOFQuantizationBinSize=39.0625;
            else
                TOFQuantizationBinSize=Gopt.TOFQuantizationBinSize;
            end
            if(~isfield(Gopt,'TOFrange'))
                TOFrange=9960.9;
            else
                TOFrange=Gopt.TOFrange;
            end
            % Form input listmode data as Castor event format
            % one event (1 time field+ 1 TOF + 1 data value + 2 crystal IDs)
            eventsize=5;
            list_in=-1*ones(length(Gopt.timefield)*eventsize,1);
            list_in(1:eventsize:end)=Gopt.timefield;
            list_in(2:eventsize:end)=Gopt.TOF;
            list_in(4:eventsize:end)=Gopt.crystalID1;
            list_in(5:eventsize:end)=Gopt.crystalID2;

            pos_val=find(list_in==-1);
            list_in(pos_val)=0;
            %
            x=reshape(x,nb_vox_X,nb_vox_Y,nb_vox_Z);
            list_out=Castor_projector_OMP(single(x),single(list_in),geo_file,projdirection,...
                               nb_vox_X,nb_vox_Y,nb_vox_Z,offset_X,offset_Y,offset_Z,...
                               fovSize_X,fovSize_Y,fovSize_Z,Gopt.num_event,dataMode,TOFinputFlag,TOFResolutionInPs,...
                               TOFQuantizationBinSize,TOFrange);
            y=list_out(pos_val);
            y=y(:);
        end
        
end
