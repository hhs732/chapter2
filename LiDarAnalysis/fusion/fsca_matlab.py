clear
clc
clf

[Z_all,R] = arcgridread('fSCA_test__all_returns_all_metrics_elevation_total_count.asc');
[Z_2m,R] = arcgridread('fSCA_test__all_returns_all_metrics_elevation_all_cover_count.asc');
[Z_nonground,R] = arcgridread('fSCA_test__all_returns_all_metrics_elevation_count.asc');

%%sort returns
Z_ground=Z_all-Z_nonground;
Z_tree=Z_2m;
Z_lowbranch=Z_nonground-Z_2m;

%create classification raster
class_raster=Z_2m*NaN;
hold_it=find(Z_tree>0 & Z_lowbranch>0);
class_raster(hold_it)=1; %Canopy with low branches
hold_it=find(Z_tree>0 & Z_lowbranch==0);
class_raster(hold_it)=2; %Canopy No low branches
hold_it=find(Z_tree==0 & Z_lowbranch>0);
class_raster(hold_it)=3; %Open with low branches
hold_it=find(Z_tree==0 & Z_lowbranch==0 & Z_ground>0);
class_raster(hold_it)=4; %Open with low branches


%% plot
subplot(1,3,1)
imagesc(Z_ground)
caxis([0 30])
subplot(1,3,2)
imagesc(Z_tree)
caxis([0 30])
subplot(1,3,3)
imagesc(Z_lowbranch)
caxis([0 30])

clf
imagesc(class_raster)
caxis([0 4])
colormap(jet(4))