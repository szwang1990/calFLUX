dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_FIO-ESM-2-0_historical_r1i1p1f1_gn.nc';
ncdisp([dname,fname]);
area = ncread([dname,fname],'areacello');
lon  = ncread([dname,fname],'lon');
lat  = ncread([dname,fname],'lat');
nlon  = ncread([dname,fname],'nlon');
nlat  = ncread([dname,fname],'nlat');
lat_bnds = ncread([dname,fname],'lat_bnds');
lon_bnds = ncread([dname,fname],'lon_bnds');

lat_bnds = reshape(lat_bnds,4,[]);
lon_bnds = reshape(lon_bnds,4,[]);

%%


figure;
% patch(lon_bnds,lat_bnds,area(:));
% patch(lon_bnds,lat_bnds,area(:),'facecolor','none');
scatter(lon_bnds(:),lat_bnds(:));

axis tight;


%%
dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc';
ncdisp([dname,fname]);
dmask = '/media/swang/CMIP6_4/UKESM1-0-LL/';
fmask = 'so_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_210001-210012_2100.nc';
ncdisp([dmask,fmask]);
mask_3d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf Inf 1],[1,1,1,1]);
mask_2d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf 1   1],[1,1,1,1]);
mask_3d(~isnan(mask_3d)) = 1;
mask_3d(isnan(mask_3d))  = 0;
mask_2d(~isnan(mask_2d)) = 1;
mask_2d(isnan(mask_2d))  = 0;

area = ncread([dname,fname],'areacello');
lon  = ncread([dname,fname],'longitude');
lat  = ncread([dname,fname],'latitude');
nlon  = ncread([dname,fname],'i');
nlat  = ncread([dname,fname],'j');
lat_bnds = ncread([dname,fname],'vertices_latitude');
lon_bnds = ncread([dname,fname],'vertices_longitude');

lat_bnds = reshape(lat_bnds,4,[]);
lon_bnds = reshape(lon_bnds,4,[]);

sec_name = {'Fram_strait','Davis_strait','Bering_strait','Denmark_strait','Barents_opening'};


%%
figure;
% patch(lon_bnds,lat_bnds,area(:));
patch(lon_bnds,lat_bnds,area(:),'facecolor','none');
% scatter(lon_bnds(:),lat_bnds(:));

axis tight;


%%
figure;
m_proj('Azimuthal Equal-area','lon',0,'lat',90,'radius',35);
[x y] = m_ll2xy(lon_bnds,lat_bnds);
patch(x,y,area(:),'facecolor','none');

axis tight;
m_grid('xtick',[-180:45:180],'xticklabel',[],'ytick',[50:10:90],'yticklabel',[]);
m_coast('patch',[.8 .8 .8],'edgecolor','k');

%%
figure;
m_proj('Azimuthal Equal-area','lon',0,'lat',90,'radius',25);
[x_bnds y_bnds] = m_ll2xy(lon_bnds,lat_bnds);
[x_cent y_cent] = m_ll2xy(lon,lat);
% patch(x,y,area(:),'facecolor','none','EdgeColor','r');
patch(x,y,mask_2d(:),'EdgeColor','r');

% color_mask = ones(length(mask_2d(:)),3);
% color_mask(find(mask_2d(:)==0),:) = 0;
% patch(x,y,color_mask,'EdgeColor','r');


axis tight;
hold on;
m_scatter(lon(:),lat(:),2,'Filled');
m_grid('xtick',[-180:45:180],'xticklabel',[],'ytick',[50:10:90],'yticklabel',[],'Color','b','LineStyle','-');
% m_coast('patch',[.8 .8 .8],'edgecolor','k');
% m_coast('patch',[.9 .9 .9],'edgecolor','k');



%%
dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc';
ncdisp([dname,fname]);
dmask = '/media/swang/CMIP6_4/UKESM1-0-LL/';
fmask = 'so_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_210001-210012_2100.nc';
fout  = 'UKESM1-0-LL_sections';
ncdisp([dmask,fmask]);
mask_3d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf Inf 1],[1,1,1,1]);
mask_2d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf 1   1],[1,1,1,1]);
mask_3d(~isnan(mask_3d)) = 1;
mask_3d(isnan(mask_3d))  = 0;
mask_2d(~isnan(mask_2d)) = 1;
mask_2d(isnan(mask_2d))  = 0;

area = ncread([dname,fname],'areacello');
lon  = ncread([dname,fname],'longitude');
lat  = ncread([dname,fname],'latitude');
nlon  = ncread([dname,fname],'i');
nlat  = ncread([dname,fname],'j');
lat_bnds = ncread([dname,fname],'vertices_latitude');
lon_bnds = ncread([dname,fname],'vertices_longitude');

lat_bnds = reshape(lat_bnds,4,[]);
lon_bnds = reshape(lon_bnds,4,[]);

sec_list = {'Fram_strait','Davis_strait','Bering_strait','Denmark_strait','Barents_opening'};
sec_list = {'Fram_strait'};
sec_list = {'Barents_opening','Denmark_strait'};
sec_list = {'Bering_strait','Denmark_strait'};

for isec=1:length(sec_list)
    sec_name = char(sec_list(isec));
    
%%
figure;
m_proj('Azimuthal Equal-area','lon',0,'lat',90,'radius',28);
[x_bnds y_bnds] = m_ll2xy(lon_bnds,lat_bnds);
[x_cent y_cent] = m_ll2xy(lon,lat);
% patch(x,y,area(:),'facecolor','none','EdgeColor','r');
patch(x_bnds,y_bnds,mask_2d(:),'EdgeColor','r');

% color_mask = ones(length(mask_2d(:)),3);
% color_mask(find(mask_2d(:)==0),:) = 0;
% patch(x,y,color_mask,'EdgeColor','r');


axis tight;
hold on;
m_scatter(lon(:),lat(:),2,'Filled');
m_grid('xtick',[-180:45:180],'xticklabel',[],'ytick',[50:10:90],'yticklabel',[],'Color','b','LineStyle','-');
% m_coast('patch',[.8 .8 .8],'edgecolor','k');
% m_coast('patch',[.9 .9 .9],'edgecolor','k');
title(['Zoom in for ', strrep(sec_name,'_','-'), ' and press ENTER']);
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
        scatter(x_cent(ind),y_cent(ind),20,'k','Filled');
        loc = [loc ind];
        clear x_sec y_sec ind;
    end
end
% scatter(x_cent(ind),y_cent(ind),30,'k','Filled');
eval([sec_name,'=[',num2str(loc),']']);
save('test.mat',sec_name,'-append');
end