clear;
clc;
close all;


h5files=dir('./*.h5')
N=numel(h5files);
startFrames=zeros(1,N);
for i=1:N
    startFrames(i)=str2double(h5files(i).name(10:end-52));
end
[~,I]=sort(startFrames,'ascend');
h5files=h5files(I);

trajectories=cell(1,N);
for i=1:N
    trajectories{i}=getTraj(h5files(i).name);
end
trajName=h5files(1).name(); % name the .mat
save([trajName,'.mat'],'trajectories','h5files');


function trajectories=getTraj(h5file)
    table=h5read(h5file,'/df_with_missing/table');
    trajData=table.values_block_0;
    [Ntraj,~]=size(trajData);
    Ntraj=floor(Ntraj/3);
    trajectories=[];
    for i=1:Ntraj
        trajectory.x=trajData((i-1)*3+1,:);
        trajectory.y=trajData((i-1)*3+2,:);
        trajectory.prob=trajData((i-1)*3+3,:);
        trajectories=[trajectories;trajectory];
    end
end