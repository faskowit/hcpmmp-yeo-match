
clc
clearvars

%% add stuff

projdir = pwd ;
addpath(genpath(projdir))

%% add read cifti
 
ciftipath = '/home/spornslab/joshstuff/git_pull/cifti-matlab' ;
addpath(ciftipath)

%% load parcs

hcpmmp_file_lh = './data/Glasser_et_al_2016_HCP_MMP1.0/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.dlabel.nii' ;
hcpmmp_file_rh = './data/Glasser_et_al_2016_HCP_MMP1.0/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.dlabel.nii' ;

hcpmmp = struct() ;
hcpmmp.lh = cifti_read(hcpmmp_file_lh) ;
hcpmmp.rh = cifti_read(hcpmmp_file_rh) ;
hcpmmp.lh.inds = hcpmmp.lh.diminfo{1}.models{1}.vertlist +1 ;
hcpmmp.rh.inds = hcpmmp.rh.diminfo{1}.models{1}.vertlist +1 ;

%%

sch_scales = 100:100:1000 ;

lh_yeo = nan(180,length(sch_scales)) ;
rh_yeo = nan(180,length(sch_scales)) ;

%%

for lll = 1:length(sch_scales) 

    sch_file = [ './data/HCP/fslr32k/cifti/Schaefer2018_' ...
        num2str(sch_scales(lll)) 'Parcels_17Networks_order.dlabel.nii' ] ;
    disp(sch_file)
    schdata = cifti_read(sch_file) ;
    lh_ind = schdata.diminfo{1}.models{1}.start:(schdata.diminfo{1}.models{1}.start+schdata.diminfo{1}.models{1}.count-1) ;
    rh_ind = schdata.diminfo{1}.models{2}.start:(schdata.diminfo{1}.models{2}.start+schdata.diminfo{1}.models{2}.count-1) ;

    lh_yeolabs = schdata.cdata(lh_ind) ;
    rh_yeolabs = schdata.cdata(rh_ind) ; 

    lh_yeolab_redu = lh_yeolabs(hcpmmp.lh.inds) ;
    rh_yeolab_redu = rh_yeolabs(hcpmmp.rh.inds) ;

    %% convert schaefer keys into network names

    sch_reg = [ schdata.diminfo{2}.maps.table.key ] ;
    sch_names = { schdata.diminfo{2}.maps.table.name } ;
    sch_yeonet = nan(length(sch_reg),1) ;

    yeonames = { '???' 'VisCent' 'VisPeri' 'SomMotA' 'SomMotB' ...
        'DorsAttnA' 'DorsAttnB' 'SalVentAttnA' 'SalVentAttnB' 'LimbicA'...
        'LimbicB' 'ContA' 'ContB' 'ContC' 'DefaultA' 'DefaultB' 'DefaultC' ...
        'TempPar' } ;

    for iii = 1:length(yeonames)

        currlab = yeonames{iii} ;

        sch_yeonet(~cellfun(@isempty,regexp(sch_names,currlab))) = iii ; 

    end

    %% do the matching

    thr_o = 80 ; 

    for iii = 1:180
        %% lh

        % get the schaefer regions here
        s_reg = lh_yeolab_redu(hcpmmp.lh.cdata==iii) ; 

        tab = tabulate(s_reg) ;
        [~,sss] = sort(tab(:,3),'descend') ;
        tab = tab(sss,:) ;

        thrind = find(cumsum(tab(:,3))>thr_o,1) ;
        s_reg_ovrlp = tab(1:thrind,1) ;

        mmo = mode(sch_yeonet(s_reg_ovrlp+1)) ; % add one to make it ind
        lh_yeo(iii,lll) = mmo(1) ;

        %% rh

            % get the schaefer regions here
        s_reg = rh_yeolab_redu(hcpmmp.rh.cdata==iii) ; 

        tab = tabulate(s_reg) ;
        [~,sss] = sort(tab(:,3),'descend') ;
        tab = tab(sss,:) ;

        thrind = find(cumsum(tab(:,3))>thr_o,1) ;
        s_reg_ovrlp = tab(1:thrind,1) ;

    %     n_ovrlp_reg = length(s_reg_ovrlp) ;
        mmo = mode(sch_yeonet(s_reg_ovrlp+1)) ; % add one to make it ind
        rh_yeo(iii,lll) = mmo(1) ;

    end

end

%% and finally a mode

% yeonames = { '???' 'VisCent' 'VisPeri' 'SomMotA' 'SomMotB' ...
%     'DorsAttnA' 'DorsAttnB' 'SalVentAttnA' 'SalVentAttnB' 'LimbicA'...
%     'LimbicB' 'ContA' 'ContB' 'ContC' 'DefaultA' 'DefaultB' 'DefaultC' ...
%     'TempPar' } ;
lh_mult_hcpyeo = mode(lh_yeo,2) ;
rh_mult_hcpyeo = mode(rh_yeo,2) ;

outtable = table() ;
nl = { hcpmmp.lh.diminfo{2}.maps.table.name } ;
nr = { hcpmmp.rh.diminfo{2}.maps.table.name } ;
outtable.names = [ nl(2:end)' ; nr(2:end)' ]; 
outtable.yeo17 = [ yeonames(lh_mult_hcpyeo)' ; yeonames(rh_mult_hcpyeo)' ] ;
outtable.yeo17ind = [ lh_mult_hcpyeo(:) ; rh_mult_hcpyeo(:) ] ;

writetable(outtable,'./hcpyeo.csv')