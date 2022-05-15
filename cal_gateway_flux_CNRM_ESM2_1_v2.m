S_ref = 34.8;
T_ref = 0;

dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_CNRM-ESM2-1_historical_r1i1p1f2_gn.nc';
lat_bnds = ncread([dname,fname],'bounds_lat');
lon_bnds = ncread([dname,fname],'bounds_lon');

dir_root = '/media/swang/CMIP6_2/';
mdl_name = 'CNRM-ESM2-1';
dir_mdl = [dir_root,mdl_name,'/'];

fsec = '~/Data/grid_and_mask/areacello_cmip6/CNRM-ESM2-1_sections.mat';
load(fsec);
scenario_list = {'historical','ssp126','ssp245','ssp370','ssp585'};
scenario_list = {'historical','ssp245','ssp585'};

ensemble_list = {'r1i1p1f2'};
year_s = 1950;
% year_e = 1948;
year_e = 2100;
rEarth = 12735000/2.0;% m

for iensem=1:length(ensemble_list)
    for iscenario=1:length(scenario_list)
         
        for iyear=year_s:year_e
            if (iyear<=2014 & ~strcmp(char(scenario_list(iscenario)),'historical')) continue;end;
            if (iyear>2014  & strcmp(char(scenario_list(iscenario)),'historical')) continue;end;
            disp('___Calculating Section tranport__________________');
            disp(['___',char(scenario_list(iscenario)),' : ',num2str(iyear)]);
            fu = ['uo_Omon_',mdl_name,'_',char(scenario_list(iscenario)),'_',char(ensemble_list(iensem)),'_*_',num2str(iyear),'.nc'];
            fv = ['vo_Omon_',mdl_name,'_',char(scenario_list(iscenario)),'_',char(ensemble_list(iensem)),'_*_',num2str(iyear),'.nc'];
            fs = ['so_Omon_',mdl_name,'_',char(scenario_list(iscenario)),'_',char(ensemble_list(iensem)),'_*_',num2str(iyear),'.nc'];
            ft = ['thetao_Omon_',mdl_name,'_',char(scenario_list(iscenario)),'_',char(ensemble_list(iensem)),'_*_',num2str(iyear),'.nc'];
            aa = dir([dir_mdl,fu]);
            fu = [dir_mdl,'/',aa.name];
            aa = dir([dir_mdl,fv]);
            fv = [dir_mdl,'/',aa.name];
            aa = dir([dir_mdl,fs]);
            fs = [dir_mdl,'/',aa.name];
            aa = dir([dir_mdl,ft]);
            ft = [dir_mdl,'/',aa.name];
            
            dep = ncread(fu,'lev');
            dep_bnds = ncread(fu,'lev_bounds');
            dz = dep_bnds(2,:) - dep_bnds(1,:);
            u3d = ncread(fu,'uo');
            v3d = ncread(fv,'vo');
            s3d = ncread(fs,'so');
            t3d = ncread(ft,'thetao');
            u3d = mean(u3d,4);
            v3d = mean(v3d,4);
            s3d = mean(s3d,4);
            t3d = mean(t3d,4);
            
            u3d = reshape(u3d,[],length(dep));
            v3d = reshape(v3d,[],length(dep));
            s3d = reshape(s3d,[],length(dep));
            t3d = reshape(t3d,[],length(dep));
           %% Fram Strait
            u2d = u3d(Fram_strait,:);
            v2d = v3d(Fram_strait,:);
            s2d = s3d(Fram_strait,:);
            t2d = t3d(Fram_strait,:);
            vol_fram = zeros(3,1);
            fw_fram  = zeros(3,1);
            ht_fram  = zeros(3,1);
            s_fram   = 0;
            t_fram   = 0;

            dx = [];
            for inod = 1:length(Fram_strait)
                nod_lons = lon_bnds(:,Fram_strait(inod));
                nod_lats = lat_bnds(:,Fram_strait(inod));
                lon1 = nod_lons(1); lat1 = nod_lats(1);
                lon2 = nod_lons(2); lat2 = nod_lats(2);
                dx = [dx rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
            end
            
            dx = dx';
            dxx = repmat(dx,1,length(dep));
            dzz = repmat(dz,length(Fram_strait),1);
            dA  = dxx.*dzz;
            loc1 = find(v2d>0);
            loc2 = find(v2d<0);
            
%             s_fram(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA));% net
%             s_fram(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             s_fram(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             t_fram(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA));% net
%             t_fram(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             t_fram(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             s_fram = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
%             t_fram = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            s_fram(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            s_fram(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            s_fram(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
            t_fram(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            t_fram(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            t_fram(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
            vol_fram(1,1) = nansum(nansum(v2d.*dA));% net
            vol_fram(2,1) = nansum(nansum(v2d(loc1).*dA(loc1)));% in
            vol_fram(3,1) = nansum(nansum(v2d(loc2).*dA(loc2)));% out
            fw_fram(1,1)  = nansum(nansum(v2d.*dA.*(S_ref-s2d)/S_ref));% net
            fw_fram(2,1)  = nansum(nansum(v2d(loc1).*dA(loc1).*(S_ref-s2d(loc1))/S_ref));% 
            fw_fram(3,1)  = nansum(nansum(v2d(loc2).*dA(loc2).*(S_ref-s2d(loc2))/S_ref));%
            ht_fram(1,1)  = 4.2e6*nansum(nansum(v2d.*dA.*(t2d-T_ref)));%
            ht_fram(2,1)  = 4.2e6*nansum(nansum(v2d(loc1).*dA(loc1).*(t2d(loc1)-T_ref)));%
            ht_fram(3,1)  = 4.2e6*nansum(nansum(v2d(loc2).*dA(loc2).*(t2d(loc2)-T_ref)));%
            
            %% Bering Strait
            u2d = u3d(Bering_strait,:);
            v2d = v3d(Bering_strait,:);
            s2d = s3d(Bering_strait,:);
            t2d = t3d(Bering_strait,:);
            vol_bering = zeros(3,1);
            fw_bering  = zeros(3,1);
            ht_bering  = zeros(3,1);
            s_bering   = 0;
            t_bering   = 0;
            dx = [];
            for inod = 1:length(Bering_strait)
                nod_lons = lon_bnds(:,Bering_strait(inod));
                nod_lats = lat_bnds(:,Bering_strait(inod));
                lon1 = nod_lons(3); lat1 = nod_lats(3);
                lon2 = nod_lons(4); lat2 = nod_lats(4);
                dx = [dx rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
            end
            
            dx = dx';
            dxx = repmat(dx,1,length(dep));
            dzz = repmat(dz,length(Bering_strait),1);
            dA  = dxx.*dzz;
            loc1 = find(v2d>0);
            loc2 = find(v2d<0);
%             s_bering(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA));% net
%             s_bering(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             s_bering(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             t_bering(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA));% net
%             t_bering(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             t_bering(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             s_bering = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
%             t_bering = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net

            s_bering(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            s_bering(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA((loc1))));% in
            s_bering(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA((loc2))));% in
            t_bering(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            t_bering(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            t_bering(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA((loc2))));% out
            vol_bering(1,1) = nansum(nansum(v2d.*dA));% net
            vol_bering(2,1) = nansum(nansum(v2d(loc1).*dA(loc1)));% in
            vol_bering(3,1) = nansum(nansum(v2d(loc2).*dA(loc2)));% out
            fw_bering(1,1)  = nansum(nansum(v2d.*dA.*(S_ref-s2d)/S_ref));%
            fw_bering(2,1)  = nansum(nansum(v2d(loc1).*dA(loc1).*(S_ref-s2d(loc1))/S_ref));%
            fw_bering(3,1)  = nansum(nansum(v2d(loc2).*dA(loc2).*(S_ref-s2d(loc2))/S_ref));%
            ht_bering(1,1)  = 4.2e6*nansum(nansum(v2d.*dA.*(t2d-T_ref)));%
            ht_bering(2,1)  = 4.2e6*nansum(nansum(v2d(loc1).*dA(loc1).*(t2d(loc1)-T_ref)));%
            ht_bering(3,1)  = 4.2e6*nansum(nansum(v2d(loc2).*dA(loc2).*(t2d(loc2)-T_ref)));%
             
            %% Davis Strait
            u2d = u3d(Davis_strait,:);
            v2d = v3d(Davis_strait,:);
            s2d = s3d(Davis_strait,:);
            t2d = t3d(Davis_strait,:);
            vol_davis = zeros(3,1);
            fw_davis  = zeros(3,1);
            ht_davis  = zeros(3,1);
            s_davis   = 0;
            t_davis   = 0;
            
            dx = [];
            for inod = 1:length(Davis_strait)
                nod_lons = lon_bnds(:,Davis_strait(inod));
                nod_lats = lat_bnds(:,Davis_strait(inod));
                lon1 = nod_lons(3); lat1 = nod_lats(3);
                lon2 = nod_lons(4); lat2 = nod_lats(4);
                dx = [dx rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
            end
            
            dx = dx';
            dxx = repmat(dx,1,length(dep));
            dzz = repmat(dz,length(Davis_strait),1);
            dA  = dxx.*dzz;
            loc1 = find(v2d>0);
            loc2 = find(v2d<0);
%             s_davis(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA));% net
%             s_davis(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             s_davis(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             t_davis(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA));% net
%             t_davis(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
%             t_davis(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             s_davis = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
%             t_davis = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            s_davis(1,1) = nansum(nansum(s2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            s_davis(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            s_davis(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
            t_davis(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA(~isnan(s2d))));% net
            t_davis(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            t_davis(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
            vol_davis(1,1) = nansum(nansum(v2d.*dA));% net
            vol_davis(2,1) = nansum(nansum(v2d(loc1).*dA(loc1)));% in
            vol_davis(3,1) = nansum(nansum(v2d(loc2).*dA(loc2)));% out
            fw_davis(1,1)  = nansum(nansum(v2d.*dA.*(S_ref-s2d)/S_ref));%
            fw_davis(2,1)  = nansum(nansum(v2d(loc1).*dA(loc1).*(S_ref-s2d(loc1))/S_ref));%
            fw_davis(3,1)  = nansum(nansum(v2d(loc2).*dA(loc2).*(S_ref-s2d(loc2))/S_ref));%
            ht_davis(1,1)  = 4.2e6*nansum(nansum(v2d.*dA.*(t2d-T_ref)));%
            ht_davis(2,1)  = 4.2e6*nansum(nansum(v2d(loc1).*dA(loc1).*(t2d(loc1)-T_ref)));%
            ht_davis(3,1)  = 4.2e6*nansum(nansum(v2d(loc2).*dA(loc2).*(t2d(loc2)-T_ref)));%
                      
            %% Denmark Strait
            u2d = u3d(Denmark_strait,:);
            v2d = v3d(Denmark_strait,:);
            s2d = s3d(Denmark_strait,:);
            t2d = t3d(Denmark_strait,:);
            vol_denmark = zeros(3,1);
            fw_denmark  = zeros(3,1);
            ht_denmark  = zeros(3,1);
            s_denmark   = 0;
            t_denmark   = 0;

            dx = []; dy = [];
            for inod = 1:length(Denmark_strait)
                nod_lons = lon_bnds(:,Denmark_strait(inod));
                nod_lats = lat_bnds(:,Denmark_strait(inod));
                lon1 = nod_lons(2); lat1 = nod_lats(2);
                lon2 = nod_lons(1); lat2 = nod_lats(1);
                dx = [dx rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
                 
                lon1 = nod_lons(2); lat1 = nod_lats(2);
                lon2 = nod_lons(3); lat2 = nod_lats(3);
                dy = [dy rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
            end
            
            % dy = flipud(dy');
            dx = dx'; dy = dy';
            dxx = repmat(dx,1,length(dep));
            dyy = repmat(dy,1,length(dep));
            dzz = repmat(dz,length(Denmark_strait),1);
            dAx  = dxx.*dzz;
            dAy  = dyy.*dzz;
            loc1u = find(u2d>0);
            loc2u = find(u2d<0);
            loc1v = find(v2d>0);
            loc2v = find(v2d<0);
            
%             s_denmark = nansum(nansum(s2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
%             t_denmark = nansum(nansum(t2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
            s_denmark(1,1) = nansum(nansum(s2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
            s_denmark(2,1) = (nansum(nansum(s2d(loc1u).*dAy(loc1u)))+nansum(nansum(s2d(loc1v).*dAx(loc1v))))/(nansum(nansum(dAy(loc1u)))+nansum(nansum(dAx(loc1v))));% in
            s_denmark(3,1) = (nansum(nansum(s2d(loc2u).*dAy(loc2u)))+nansum(nansum(s2d(loc2v).*dAx(loc2v))))/(nansum(nansum(dAy(loc2u)))+nansum(nansum(dAx(loc2v))));% out
            t_denmark(1,1) = nansum(nansum(t2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(t2d))+dAy(~isnan(s2d))));% net
            t_denmark(2,1) = (nansum(nansum(t2d(loc1u).*dAy(loc1u)))+nansum(nansum(t2d(loc1v).*dAx(loc1v))))/(nansum(nansum(dAy(loc1u)))+nansum(nansum(dAx(loc1v))));% in
            t_denmark(3,1) = (nansum(nansum(t2d(loc2u).*dAy(loc2u)))+nansum(nansum(t2d(loc2v).*dAx(loc2v))))/(nansum(nansum(dAy(loc2u)))+nansum(nansum(dAx(loc2v))));% out
            vol_denmark(1,1) = nansum(nansum(u2d.*dAy+v2d.*dAx));
            vol_denmark(2,1) = nansum(nansum(u2d(loc1u).*dAy(loc1u)))+nansum(nansum(v2d(loc1v).*dAx(loc1v)));
            vol_denmark(3,1) = nansum(nansum(u2d(loc2u).*dAy(loc2u)))+nansum(nansum(v2d(loc2v).*dAx(loc2v)));
            fw_denmark(1,1)  = nansum(nansum(u2d.*dAy.*(S_ref-s2d)/S_ref))+nansum(nansum(v2d.*dAx.*(S_ref-s2d)/S_ref));%
            fw_denmark(2,1)  = nansum(nansum(u2d(loc1u).*dAy(loc1u).*(S_ref-s2d(loc1u))/S_ref))+nansum(nansum(v2d(loc1v).*dAx(loc1v).*(S_ref-s2d(loc1v))/S_ref));%
            fw_denmark(3,1)  = nansum(nansum(u2d(loc2u).*dAy(loc2u).*(S_ref-s2d(loc2u))/S_ref))+nansum(nansum(v2d(loc2v).*dAx(loc2v).*(S_ref-s2d(loc2v))/S_ref));%
            ht_denmark(1,1)  = 4.2e6*nansum(nansum(u2d.*dAy.*(t2d-T_ref)))+nansum(nansum(v2d.*dAx.*(t2d-T_ref)));%
            ht_denmark(2,1)  = 4.2e6*nansum(nansum(u2d(loc1u).*dAy(loc1u).*(t2d(loc1u)-T_ref)))+nansum(nansum(v2d(loc1v).*dAx(loc1v).*(t2d(loc1v)-T_ref)));%
            ht_denmark(3,1)  = 4.2e6*nansum(nansum(u2d(loc2u).*dAy(loc2u).*(t2d(loc2u)-T_ref)))+nansum(nansum(v2d(loc2v).*dAx(loc2v).*(t2d(loc2v)-T_ref)));%

            % voltransp_denmark = nansum(nansum(u2d.*dA));%
            
            %% Barents Sea Opening
            u2d = u3d(Barents_opening,:);
            v2d = v3d(Barents_opening,:);
            s2d = s3d(Barents_opening,:);
            t2d = t3d(Barents_opening,:);
            vol_bso = zeros(3,1);
            fw_bso  = zeros(3,1);
            ht_bso  = zeros(3,1);
            s_bso   = 0;
            t_bso   = 0;
            dx = []; dy = [];
            for inod = 1:length(Barents_opening)
                nod_lons = lon_bnds(:,Barents_opening(inod));
                nod_lats = lat_bnds(:,Barents_opening(inod));
                lon1 = nod_lons(1); lat1 = nod_lats(1);
                lon2 = nod_lons(2); lat2 = nod_lats(2);                
                dx = [dx rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
                
                lon1 = nod_lons(3); lat1 = nod_lats(3);
                lon2 = nod_lons(2); lat2 = nod_lats(2);
                dy = [dy rEarth*acos(cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)+sind(lat1)*sind(lat2))];
                
            end
            
            dx = dx'; dy = dy';
            dxx = repmat(dx,1,length(dep));
            dyy = repmat(dy,1,length(dep));
            dzz = repmat(dz,length(Barents_opening),1);
            dAx  = dxx.*dzz;
            dAy  = dyy.*dzz;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!
            % the calculation of BSO is grid dependent
            % In FIO-ESM, it is BSO inflow when uo>0 or vo<0
            loc1u = find(u2d>0);
            loc2u = find(u2d<0);
            loc1v = find(v2d>0);
            loc2v = find(v2d<0);

            % s_bso(1,1) = nansum(nansum(s2d.*(dAx+dAy)))/nansum(nansum(dAx+dAy));% net
            % s_bso(2,1) = nansum(nansum(s2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            % s_bso(3,1) = nansum(nansum(s2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
            % t_bso(1,1) = nansum(nansum(t2d.*dA))/nansum(nansum(dA));% net
            % t_bso(2,1) = nansum(nansum(t2d(loc1).*dA(loc1)))/nansum(nansum(dA(loc1)));% in
            % t_bso(3,1) = nansum(nansum(t2d(loc2).*dA(loc2)))/nansum(nansum(dA(loc2)));% out
%             s_bso = nansum(nansum(s2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
%             t_bso = nansum(nansum(t2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net

            s_bso(1,1) = nansum(nansum(s2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
            s_bso(2,1) = (nansum(nansum(s2d(loc1u).*dAy(loc1u)))+nansum(nansum(s2d(loc1v).*dAx(loc1v))))/(nansum(nansum(dAy(loc1u)))+nansum(nansum(dAx(loc1v))));% in
            s_bso(3,1) = (nansum(nansum(s2d(loc2u).*dAy(loc2u)))+nansum(nansum(s2d(loc2v).*dAx(loc2v))))/(nansum(nansum(dAy(loc2u)))+nansum(nansum(dAx(loc2v))));% out
            t_bso(1,1) = nansum(nansum(t2d.*(dAx+dAy)))/nansum(nansum(dAx(~isnan(s2d))+dAy(~isnan(s2d))));% net
            t_bso(2,1) = (nansum(nansum(t2d(loc1u).*dAy(loc1u)))+nansum(nansum(t2d(loc1v).*dAx(loc1v))))/(nansum(nansum(dAy(loc1u)))+nansum(nansum(dAx(loc1v))));% in
            t_bso(3,1) = (nansum(nansum(t2d(loc2u).*dAy(loc2u)))+nansum(nansum(t2d(loc2v).*dAx(loc2v))))/(nansum(nansum(dAy(loc2u)))+nansum(nansum(dAx(loc2v))));% in
            vol_bso(1,1) = nansum(nansum(u2d.*dAy+v2d.*dAx));
            vol_bso(2,1) = nansum(nansum(u2d(loc1u).*dAy(loc1u)))+nansum(nansum(v2d(loc1v).*dAx(loc1v)));
            vol_bso(3,1) = nansum(nansum(u2d(loc2u).*dAy(loc2u)))+nansum(nansum(v2d(loc2v).*dAx(loc2v)));
            fw_bso(1,1)  = nansum(nansum(u2d.*dAy.*(S_ref-s2d)/S_ref))+nansum(nansum(v2d.*dAx.*(S_ref-s2d)/S_ref));%
            fw_bso(2,1)  = nansum(nansum(u2d(loc1u).*dAy(loc1u).*(S_ref-s2d(loc1u))/S_ref))+nansum(nansum(v2d(loc1v).*dAx(loc1v).*(S_ref-s2d(loc1v))/S_ref));%
            fw_bso(3,1)  = nansum(nansum(u2d(loc2u).*dAy(loc2u).*(S_ref-s2d(loc2u))/S_ref))+nansum(nansum(v2d(loc2v).*dAx(loc2v).*(S_ref-s2d(loc2v))/S_ref));%
            ht_bso(1,1)  = 4.2e6*nansum(nansum(u2d.*dAy.*(t2d-T_ref)))+nansum(nansum(v2d.*dAx.*(t2d-T_ref)));%
            ht_bso(2,1)  = 4.2e6*nansum(nansum(u2d(loc1u).*dAy(loc1u).*(t2d(loc1u)-T_ref)))+nansum(nansum(v2d(loc1v).*dAx(loc1v).*(t2d(loc1v)-T_ref)));%
            ht_bso(3,1)  = 4.2e6*nansum(nansum(u2d(loc2u).*dAy(loc2u).*(t2d(loc2u)-T_ref)))+nansum(nansum(v2d(loc2v).*dAx(loc2v).*(t2d(loc2v)-T_ref)));%
            
            % voltransp_bso = nansum(nansum(u2d.*dAy-v2d.*dAx));%
            
            fout = [dir_mdl,'/Gateway_flux_yearly_',char(scenario_list(iscenario)),'_',char(ensemble_list(iensem)),'_',num2str(iyear)];
            save(fout,'vol_fram','fw_fram','ht_fram','s_fram','t_fram',...
                      'vol_bering','fw_bering','ht_bering','s_bering','t_bering',...
                      'vol_davis','fw_davis','ht_davis','s_davis','t_davis',...
                      'vol_denmark','fw_denmark','ht_denmark','s_denmark','t_denmark',...
                      'vol_bso','fw_bso','ht_bso','s_bso','t_bso');
        end
    end
end

            
            %%
            return;
            
            scatter(x_cent(Fram_strait),y_cent(Fram_strait),20,'k','Filled');
            m_scatter(lon_bnds(Fram_strait(1)),lat_bnds(Fram_strait(1)),20,'r','Filled');
            
            figure;
            contourf(dA);
            colorbar;
            
            figure;
            contourf(s2d);
            colorbar;
            
            figure;
            contourf(u2d);
            colorbar;
            colormap(cmocean('balance',10));
            caxis([-0.1 0.1]);
            
            figure;
            contourf(v2d);
            colorbar;
            colormap(cmocean('balance',10));
            caxis([-0.1 0.1]);
            
