clear;
clc;
close all;


h5files=dir('./*.h5');
N=numel(h5files);
startFrames=zeros(1,N);
for i=1:N
    startFrames(i)=str2double(h5files(i).name(10:end-52));
end
[~,I]=sort(startFrames,'ascend');
h5files=h5files(I);
RJONOFF=[0,1,0,1,1,0,1,1,0];
enable= [1,1,1,0,1,1,0,1,0];
th=0.95;

figure;
hold on;
distOn=0;
for i=1:N
    if RJONOFF(i)&&enable(i)
        distOn=distOn+DispTraj(h5files(i).name,th,'b');
    end
end
fprintf('RJON: traveled %.2f\n',distOn);

figure;
hold on;
distOff=0;
for i=1:N
    if ~RJONOFF(i)&&enable(i)
        distOff=distOff+DispTraj(h5files(i).name,th,'b');
    end
end
fprintf('RJOFF: traveled %.4f\n',distOff);
fprintf('RJOFF: traveled %.4f m\n',distOff/500*0.25); % distance [m]

function dist=DispTraj(h5file,th,color)
    table=h5read(h5file,'/df_with_missing/table');
    trajData=table.values_block_0;
    prob=trajData(3,:);
    x=trajData(1,:);
    y=trajData(2,:);
    delta=sqrt(diff(x).^2+diff(y).^2);
    filter=delta>mean(delta)+3*std(delta);
    filter=logical([0,filter]);
    x(prob<th)=[];
    y(prob<th)=[];
    dx=diff(x);
    dy=diff(y);
    traj=x+1i*y;
    plot(traj,color);
    dist=sum(sqrt(dx.^2+dy.^2));
end