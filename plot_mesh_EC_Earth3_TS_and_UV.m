%%
mdl_name = 'EC-Earth3';
dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_EC-Earth3_historical_r1i1p1f1_gn.nc';
fout  = [mdl_name,'_sections_v1.mat'];
% ncdisp([dname,fname]);
dmask = '/media/swang/CMIP6_2/EC-Earth3/';
fmask = 'so_Omon_EC-Earth3_ssp585_r1i1p1f1_gn_210001-210012_2100.nc';
% ncdisp([dmask,fmask]);
mask_3d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf Inf 1],[1,1,1,1]);
mask_2d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf 1   1],[1,1,1,1]);
mask_3d(~isnan(mask_3d)) = 1;
mask_3d(isnan(mask_3d))  = 0;
mask_2d(~isnan(mask_2d)) = 1;
mask_2d(isnan(mask_2d))  = 0;

area = ncread([dname,fname],'areacello');
lon  = ncread([dname,fname],'longitude');
lat  = ncread([dname,fname],'latitude');
nlon = ncread([dname,fname],'i');
nlat = ncread([dname,fname],'j');
fname = 'EC-Earth3-nemo-vertices-ORCA1-t-grid.nc';
lat_bnds = ncread([dname,fname],'vertices_latitude');
lon_bnds = ncread([dname,fname],'vertices_longitude');
% ==========================================================
lat_bnds = reshape(lat_bnds,4,[]);
lon_bnds = reshape(lon_bnds,4,[]);

% ======================
fu   = 'uo_Omon_EC-Earth3_historical_r1i1p1f1_gn_195001-195012_1950.nc';
fv   = 'vo_Omon_EC-Earth3_historical_r1i1p1f1_gn_195001-195012_1950.nc';

lon_u  = ncread([dmask,fu],'longitude');
lat_u  = ncread([dmask,fu],'latitude');
lat_bnds_u = ncread([dmask,fu],'vertices_latitude');
lon_bnds_u = ncread([dmask,fu],'vertices_longitude');
lat_bnds_u = reshape(lat_bnds_u,4,[]);
lon_bnds_u = reshape(lon_bnds_u,4,[]);

lon_v  = ncread([dmask,fv],'longitude');
lat_v  = ncread([dmask,fv],'latitude');
lat_bnds_v = ncread([dmask,fv],'vertices_latitude');
lon_bnds_v = ncread([dmask,fv],'vertices_longitude');
lat_bnds_v = reshape(lat_bnds_v,4,[]);
lon_bnds_v = reshape(lon_bnds_v,4,[]);
% ======================
sec_list = {'Fram_strait','Davis_strait','Bering_strait','Denmark_strait','Barents_opening'};
% sec_list = {'Fram_strait'};
% sec_list = {'Barents_opening','Denmark_strait'};
% sec_list = {'Bering_strait','Denmark_strait'};
fout  = [mdl_name,'_sections_TS_and_UV.mat'];
save([dname,fout],'mdl_name');

%%
for isec=1:length(sec_list)
    sec_name = char(sec_list(isec));
    

figure;
m_proj('Azimuthal Equal-area','lon',0,'lat',90,'radius',35);

% m_proj('Azimuthal Equal-area','lon',180,'lat',90,'radius',90);

[x_bnds y_bnds] = m_ll2xy(lon_bnds,lat_bnds);
[x_cent y_cent] = m_ll2xy(lon,lat);
[x_cent_u y_cent_u] = m_ll2xy(lon_u,lat_u);
[x_cent_v y_cent_v] = m_ll2xy(lon_v,lat_v);
% patch(x,y,area(:),'facecolor','none','EdgeColor','r');
patch(x_bnds,y_bnds,mask_2d(:),'EdgeColor','r');

% color_mask = ones(length(mask_2d(:)),3);
% color_mask(find(mask_2d(:)==0),:) = 0;
% patch(x,y,color_mask,'EdgeColor','r');


axis tight;
hold on;
m_scatter(lon(:),lat(:),5,'Filled','MarkerFaceColor','k');
hold on;
m_scatter(lon_u(:),lat_u(:),10,'Filled','MarkerFaceColor','b');
hold on;
m_scatter(lon_v(:),lat_v(:),10,'Filled','MarkerFaceColor','r');

m_grid('xtick',[-180:45:180],'xticklabel',[],'ytick',[-90:10:90],'yticklabel',[],'Color','b','LineStyle','-');
% m_coast('patch',[.8 .8 .8],'edgecolor','k');
% m_coast('patch',[.9 .9 .9],'edgecolor','k');
title(['Zoom in for TS grid points (BLACK) of ', strrep(sec_name,'_','-'), ' and press ENTER']);
pause;
loc = [];
while(true)
    [x_sec,y_sec] = ginput(1);
    if(isempty(x_sec))
        break;
    else
        title('Pick up grid points. Press ENTER to stop!')
        [~,ind] = min(sqrt((x_cent(:)-x_sec).^2+(y_cent(:)-y_sec).^2));
        hold on;
        scatter(x_cent(ind),y_cent(ind),30,'k','Filled');
        loc = [loc ind];
        clear x_sec y_sec ind;
    end
end

sec_out = [sec_name,'_TS'];
eval([sec_out,'=[',num2str(loc),']']);
% save('test.mat',sec_name,'-append');
save([dname,fout],sec_out,'-append');
%================================================
title(['Now select the corresponding U grid points (BLUE) for ', strrep(sec_name,'_','-'),': blue']);
loc = [];
while(true)
    [x_sec,y_sec] = ginput(1);
    if(isempty(x_sec))
        break;
    else
        title('Pick up grid points. Press ENTER to stop!')
        [~,ind] = min(sqrt((x_cent_u(:)-x_sec).^2+(y_cent_u(:)-y_sec).^2));
        hold on;
        scatter(x_cent_u(ind),y_cent_u(ind),30,'b','Filled');
        loc = [loc ind];
        clear x_sec y_sec ind;
    end
end
% scatter(x_cent(ind),y_cent(ind),30,'k','Filled');
sec_out = [sec_name,'_U'];
eval([sec_out,'=[',num2str(loc),']']);
% save('test.mat',sec_name,'-append');
save([dname,fout],sec_out,'-append');
%================================================
title(['Now select the corresponding V grid points (RED) for ', strrep(sec_name,'_','-')]);
loc = [];
while(true)
    [x_sec,y_sec] = ginput(1);
    if(isempty(x_sec))
        break;
    else
        title('Pick up grid points. Press ENTER to stop!')
        [~,ind] = min(sqrt((x_cent_v(:)-x_sec).^2+(y_cent_v(:)-y_sec).^2));
        hold on;
        scatter(x_cent_v(ind),y_cent_v(ind),30,'r','Filled');
        loc = [loc ind];
        clear x_sec y_sec ind;
    end
end
% scatter(x_cent(ind),y_cent(ind),30,'k','Filled');
sec_out = [sec_name,'_V'];
eval([sec_out,'=[',num2str(loc),']']);
% save('test.mat',sec_name,'-append');
save([dname,fout],sec_out,'-append');

end

%%
return;
so_3d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf Inf 1],[1,1,1,1]);
so_2d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf 1   1],[1,1,1,1]);
so_3d = reshape(so_3d,[],40);
bering = so_3d(Bering_strait,:);
fram = so_3d(Fram_strait,:);
barents = so_3d(Barents_opening,:);
%%
figure;
contourf(fram);
colorbar;


%% Fram Strait
hold on;
scatter(x_cent(Fram_strait),y_cent(Fram_strait),20,'b','Filled');

%
hold on;
m_scatter(lon_bnds(1,Fram_strait(1)),lat_bnds(1,Fram_strait(1)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Fram_strait(1)),lat_bnds(2,Fram_strait(1)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Fram_strait(1)),lat_bnds(3,Fram_strait(1)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Fram_strait(1)),lat_bnds(4,Fram_strait(1)),20,'m','Filled')

%
hold on;
m_scatter(lon_bnds(1,Fram_strait(3)),lat_bnds(1,Fram_strait(3)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Fram_strait(3)),lat_bnds(2,Fram_strait(3)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Fram_strait(3)),lat_bnds(3,Fram_strait(3)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Fram_strait(3)),lat_bnds(4,Fram_strait(3)),20,'m','Filled')

%% Bering Strait
hold on;
scatter(x_cent(Bering_strait),y_cent(Bering_strait),20,'k','Filled');

%
hold on;
m_scatter(lon_bnds(1,Bering_strait(1)),lat_bnds(1,Bering_strait(1)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Bering_strait(1)),lat_bnds(2,Bering_strait(1)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Bering_strait(1)),lat_bnds(3,Bering_strait(1)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Bering_strait(1)),lat_bnds(4,Bering_strait(1)),20,'m','Filled')

%% Davis Strait
hold on;
scatter(x_cent(Davis_strait),y_cent(Davis_strait),20,'k','Filled');

%
hold on;
m_scatter(lon_bnds(1,Davis_strait(1)),lat_bnds(1,Davis_strait(1)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Davis_strait(1)),lat_bnds(2,Davis_strait(1)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Davis_strait(1)),lat_bnds(3,Davis_strait(1)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Davis_strait(1)),lat_bnds(4,Davis_strait(1)),20,'m','Filled')

%% Denmark Strait
hold on;
scatter(x_cent(Denmark_strait),y_cent(Denmark_strait),20,'k','Filled');

%
hold on;
m_scatter(lon_bnds(1,Denmark_strait(1)),lat_bnds(1,Denmark_strait(1)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Denmark_strait(1)),lat_bnds(2,Denmark_strait(1)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Denmark_strait(1)),lat_bnds(3,Denmark_strait(1)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Denmark_strait(1)),lat_bnds(4,Denmark_strait(1)),20,'m','Filled')
%% Barents Sea Opening
hold on;
scatter(x_cent(Barents_opening),y_cent(Barents_opening),20,'k','Filled');

%
hold on;
m_scatter(lon_bnds(1,Barents_opening(1)),lat_bnds(1,Barents_opening(1)),20,'r','Filled')

hold on;
m_scatter(lon_bnds(2,Barents_opening(1)),lat_bnds(2,Barents_opening(1)),20,'b','Filled')

hold on;
m_scatter(lon_bnds(3,Barents_opening(1)),lat_bnds(3,Barents_opening(1)),20,'c','Filled')

hold on;
m_scatter(lon_bnds(4,Barents_opening(1)),lat_bnds(4,Barents_opening(1)),20,'m','Filled')

