%%
mdl_name = 'GISS-E2-1-G';
dname = '/home/swang/Data/grid_and_mask/areacello_cmip6/';
fname = 'areacello_Ofx_GISS-E2-1-G_piControl_r1i1p1f1_gn.nc';
fout  = [mdl_name,'_sections.mat'];
% ncdisp([dname,fname]);
dmask = '/media/swang/CMIP6_5/GISS-E2-1-G/';
fmask = 'so_Omon_GISS-E2-1-G_ssp585_r1i1p1f2_gn_208101-210012_2100.nc';
% ncdisp([dmask,fmask]);
mask_3d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf Inf 1],[1,1,1,1]);
mask_2d = ncread([dmask,fmask],'so',[1,1,1,1],[Inf Inf 1   1],[1,1,1,1]);
mask_3d(~isnan(mask_3d)) = 1;
mask_3d(isnan(mask_3d))  = 0;
mask_2d(~isnan(mask_2d)) = 1;
mask_2d(isnan(mask_2d))  = 0;

area = ncread([dname,fname],'areacello');
lon_mdl  = ncread([dname,fname],'lon');
lat_mdl  = ncread([dname,fname],'lat');
[lat lon] = meshgrid(lat_mdl,lon_mdl);

lat_bnds_orig = ncread([dname,fname],'lat_bnds');
lon_bnds_orig = ncread([dname,fname],'lon_bnds');

lat_bnds = zeros([4,size(lat)]);
lon_bnds = zeros([4,size(lat)]);
for j=1:length(lat_mdl)
    for i=1:length(lon_mdl)
        lat_bnds(1,i,j) = lat_bnds_orig(1,j);
        lat_bnds(2,i,j) = lat_bnds_orig(2,j);
        lat_bnds(3,i,j) = lat_bnds_orig(2,j);
        lat_bnds(4,i,j) = lat_bnds_orig(1,j);
        lon_bnds(1,i,j) = lon_bnds_orig(1,i);
        lon_bnds(2,i,j) = lon_bnds_orig(1,i);
        lon_bnds(3,i,j) = lon_bnds_orig(2,i);
        lon_bnds(4,i,j) = lon_bnds_orig(2,i);
    end
end

lat_bnds = reshape(lat_bnds,4,[]);
lon_bnds = reshape(lon_bnds,4,[]);

sec_list = {'Fram_strait','Davis_strait','Bering_strait','Denmark_strait','Barents_opening'};
% sec_list = {'Fram_strait'};
% sec_list = {'Barents_opening','Denmark_strait'};
% sec_list = {'Bering_strait','Denmark_strait'};
% save([dname,fout],'mdl_name');

%%
for isec=1:length(sec_list)
    sec_name = char(sec_list(isec));
    

figure;
m_proj('Azimuthal Equal-area','lon',0,'lat',90,'radius',40);
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
% save('test.mat',sec_name,'-append');
save([dname,fout],sec_name,'-append');
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
scatter(x_cent(Fram_strait),y_cent(Fram_strait),20,'k','Filled');

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

