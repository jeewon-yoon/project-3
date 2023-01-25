clear all;

% Creating .mat file list into designated directory 
filename = ['*.mat']; 
current_path = [pwd '/'];
filelist = dir([current_path filename]);
nfile = size(filelist,1);

% make FMRIraw (merging data to single file)
tmpname = [];
tmpstr = [];
tmpcell = [];
FMRIraw_Control = {};

for i = 1:nfile
    display(i);
    tmpname = filelist(i).name;
    tmpstr = load(tmpname);
    tmpcell = struct2cell(tmpstr);
    FMRIraw_Control{i} = tmpcell{1}
end

save FMRIraw_Control 


%% 
% Make SLD matrix for individuals in ADHD group
clear all;

load('FMRIraw_Control.mat')

% single linkage matrix

Original.MST = []; Original.A = []; Original.Dx = []; Original.Cx =[];

for i = 1:size(FMRIraw_Control,2)
    display(i);
    [Original.MST(:,:,i), Original.A(:,:,i), Original.Dx(:,:,i), Original.Cx(:,:,i)] = shapeofnetwork_cc_pc(FMRIraw_Control{i},'Pearson');
    SLD_Group_Control{i} = Original.Dx(:,:,i);
    Cx_Group{i} = Original. Cx(:,:,i);
    A_Group{i} = Original.A(:,:,i);
end

save Original_Control Original
save SLD_Group_Control SLD_Group_Control
save Cx_Group_Control Cx_Group
save A_Group_Control A_Group

% Draw single linkage matrix figures for 10 sample cases
for i = 21:30
    subplot(2,5,i-20); imagesc(Original.Dx(:,:,i-20));
    axis square
    set(colorbar,'Visible','on');   
    caxis([0 1]);
end
colormap(flipud('hot')); set(gcf,'Color', 'white');
saveas(gcf,'single linakge matrix, ADHD_sample10','png')
