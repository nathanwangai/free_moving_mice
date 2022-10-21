%% Neural Signal Processing Functions
classdef NSPF
    properties (Constant)
        alSize=25;
        atlSize=20;
        speedTh=5;
    end
    methods(Static)
        function Map=PlotMap(Coor,size,upsample,prjIm,colormap,contrast)
            if nargin==3
                I=zeros(size*upsample,size*upsample);
                for i=1:numel(Coor)
                    coord=Coor{i}*upsample;
                    centroid=mean(coord,2);
                    coord=coord(:);
                    I=insertShape(I,'Polygon',coord','Color','white','LineWidth',2);
                    %I=insertShape(I,'FilledPolygon',coord','Color','white','Opacity',0.5);
                end
            else
                I=imread(prjIm);
                I=imresize(I,[size*upsample,size*upsample]);
                if exist(colormap,'file')
                    low=contrast(1);
                    high=contrast(2);
                    load(colormap,'LUT');
                    I=double(I);
                    I=(I-low)/(high-low)*255;
                    I=ind2rgb(uint8(I),LUT);
                end
                for i=1:numel(Coor)
                    coord=Coor{i}*upsample;
                    centroid=mean(coord,2);
                    coord=coord(:);
                    I=insertShape(I,'Polygon',coord','Color','white','LineWidth',2);
                    %I=insertShape(I,'FilledPolygon',coord','Color','white','Opacity',0.2);
                    I=insertText(I,centroid',sprintf('%d',i),'TextColor','white',...
                        'AnchorPoint','Center','BoxOpacity',0,'FontSize',24);
                end
            end
            figure('Name','Neuron Map');
            image(I);
            axis image;
            Map=I;
        end
        function dFF_raw=GetDFF(rawTraces,BGpct)
            [~, nT]=size(rawTraces);
            F=sort(rawTraces,2,'ascend');
            F=repmat(F(:,round(nT*BGpct/100)),[1,nT]);
            dFF_raw=rawTraces./F-1;
        end
        function rawTraces=GetRawTraces(rawData,Coor)
            [nR,nC,nT]=size(rawData);
            rawDataMat=reshape(rawData,[nR*nC,nT]);
            nNeuron=numel(Coor);
            maskMat=zeros(nNeuron,nR*nC);
            for i=1:nNeuron
                poly=Coor{i};
                mask=poly2mask(poly(1,:),poly(2,:),nC,nR);
                maskMat(i,:)=reshape(mask,[1,nC*nR])/sum(sum(mask));
            end
            rawTraces=maskMat*rawDataMat;
        end
        function traj_interp=TrajInterp(traj,prob_th,diff_th,nVec)
            if nargin==3
                nVec=numel(traj.x);
            end
            filter1=traj.prob<prob_th;
            filter2=zeros(1,numel(filter1));
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                filter2(startIdx:endIdx)=([0,abs(diff(traj.x(startIdx:endIdx)))]>diff_th)...
                    |([0,abs(diff(traj.y(startIdx:endIdx)))]>diff_th);
                startIdx=startIdx+nVec(i);
            end
            filter=filter1|filter2;
            traj_filter.x=traj.x;
            traj_filter.y=traj.y;
            traj_filter.x(filter)=[];
            traj_filter.y(filter)=[];
            idx=1:numel(traj.x);
            idx_filter=idx;
            idx_filter(filter)=[];
            traj_interp.x=interp1(idx_filter,traj_filter.x,idx);
            traj_interp.y=interp1(idx_filter,traj_filter.y,idx);
        end
        function [combinedTraj,nVec]=CombineTraj(trajCell)
            traj1.x=[];
            traj1.y=[];
            traj1.prob=[];
            traj2.x=[];
            traj2.y=[];
            traj2.prob=[];
            nVec=[];
            for i=1:numel(trajCell)
                trajectory_both=trajCell{i};
                traj1.x=[traj1.x,trajectory_both(1).x];
                traj1.y=[traj1.y,trajectory_both(1).y];
                traj1.prob=[traj1.prob,trajectory_both(1).prob];
                traj2.x=[traj2.x,trajectory_both(2).x];
                traj2.y=[traj2.y,trajectory_both(2).y];
                traj2.prob=[traj2.prob,trajectory_both(2).prob];
                nVec=[nVec,numel(trajectory_both(1).x)];
            end
            combinedTraj=[traj1,traj2];
        end
        function vel=GetVelocity(combinedTraj,nVec)
            vel=zeros(1,sum(nVec));
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                vel_x=[0,diff(combinedTraj.x(startIdx:endIdx))];
                vel_y=[0,diff(combinedTraj.y(startIdx:endIdx))];
                vel(startIdx:endIdx)=sqrt(vel_x.^2+vel_y.^2);
                startIdx=startIdx+nVec(i);
            end
        end
        function angle=GetAngle(combinedTraj1,combinedTraj2,nVec)
            angle=zeros(1,sum(nVec));
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                vec_x=combinedTraj2.x(startIdx:endIdx)-combinedTraj1.x(startIdx:endIdx);
                vec_y=combinedTraj2.y(startIdx:endIdx)-combinedTraj1.y(startIdx:endIdx);
                angle(startIdx:endIdx)=angle(vec_x+1i*vec_y);
                startIdx=startIdx+nVec(i);
            end
            %angle=unwrap(angle);
        end
        function anglevel=GetAngleVel(combinedTraj1,combinedTraj2,nVec)
            anglevel=zeros(1,sum(nVec));
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                vec_x=combinedTraj2.x(startIdx:endIdx)-combinedTraj1.x(startIdx:endIdx);
                vec_y=combinedTraj2.y(startIdx:endIdx)-combinedTraj1.y(startIdx:endIdx);
                anglevel(startIdx:endIdx)=abs([diff(unwrap(angle(vec_x+1i*vec_y))),0]);
                startIdx=startIdx+nVec(i);
            end
            %angle=unwrap(angle);
        end
        function vel_ds=GetVelocity_ds(combinedTraj,nVec,sessionLen)
            vel_ds=zeros(1,numel(nVec)*sessionLen);
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                traj.x=combinedTraj.x(startIdx:endIdx);
                traj.y=combinedTraj.y(startIdx:endIdx);
                idx=(1:nVec(i))/nVec(i)*sessionLen;
                traj_ds.x=interp1(idx,traj.x,1:sessionLen);
                traj_ds.y=interp1(idx,traj.y,1:sessionLen);
                vel_ds_x=[0,diff(traj_ds.x)];
                vel_ds_y=[0,diff(traj_ds.y)];
                vel_ds((i-1)*sessionLen+1:i*sessionLen)=sqrt(vel_ds_x.^2+vel_ds_y.^2);
                startIdx=startIdx+nVec(i);
            end
        end
        function angle_ds=GetAngle_ds(combinedTraj1,combinedTraj2,nVec,sessionLen)
            angle_ds=zeros(1,numel(nVec)*sessionLen);
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                traj1.x=combinedTraj1.x(startIdx:endIdx);
                traj1.y=combinedTraj1.y(startIdx:endIdx);
                traj2.x=combinedTraj2.x(startIdx:endIdx);
                traj2.y=combinedTraj2.y(startIdx:endIdx);
                idx=(1:nVec(i))/nVec(i)*sessionLen;
                traj1_ds.x=interp1(idx,traj1.x,1:sessionLen);
                traj1_ds.y=interp1(idx,traj1.y,1:sessionLen);
                traj2_ds.x=interp1(idx,traj2.x,1:sessionLen);
                traj2_ds.y=interp1(idx,traj2.y,1:sessionLen);
                vec_ds_x=traj2_ds.x-traj1_ds.x;
                vec_ds_y=traj2_ds.y-traj1_ds.y;
                angle_ds((i-1)*sessionLen+1:i*sessionLen)=angle(vec_ds_x+1i*vec_ds_y);
                startIdx=startIdx+nVec(i);
            end
            %angle_ds=unwrap(angle_ds);
        end
        function anglevel_ds=GetAngleVel_ds(combinedTraj1,combinedTraj2,nVec,sessionLen)
            anglevel_ds=zeros(1,numel(nVec)*sessionLen);
            startIdx=1;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                traj1.x=combinedTraj1.x(startIdx:endIdx);
                traj1.y=combinedTraj1.y(startIdx:endIdx);
                traj2.x=combinedTraj2.x(startIdx:endIdx);
                traj2.y=combinedTraj2.y(startIdx:endIdx);
                idx=(1:nVec(i))/nVec(i)*sessionLen;
                traj1_ds.x=interp1(idx,traj1.x,1:sessionLen);
                traj1_ds.y=interp1(idx,traj1.y,1:sessionLen);
                traj2_ds.x=interp1(idx,traj2.x,1:sessionLen);
                traj2_ds.y=interp1(idx,traj2.y,1:sessionLen);
                vec_ds_x=traj2_ds.x-traj1_ds.x;
                vec_ds_y=traj2_ds.y-traj1_ds.y;
                anglevel_ds((i-1)*sessionLen+1:i*sessionLen)=abs([diff(unwrap(angle(vec_ds_x+1i*vec_ds_y))),0]);
                startIdx=startIdx+nVec(i);
            end
            %angle_ds=unwrap(angle_ds);
        end
        function traj_filter=GetTrajFilter(RJONOFF,nVec)
            traj_filter=[];
            for i=1:numel(RJONOFF)
                traj_filter=[traj_filter,RJONOFF(i)*ones(1,nVec(i))];
            end
        end
        function [corr,p]=PlotCorr(dFmat,valVec,title)
            [N,~]=size(dFmat);
            corr=zeros(1,N);
            p=zeros(1,N);
            for i=1:N
                [mat,pmat]=corrcoef(dFmat(i,:),valVec);
                corr(i)=mat(1,2);
                p(i)=pmat(1,2);
            end
            H=figure('Name',title);
            hold on;
            idx=1:N;
            b=bar(idx(p<0.05),corr(p<0.05));
            bmax=bar(idx(corr==max(corr)),corr(corr==max(corr)));
            b.LineWidth=1.5;
            bmax.LineWidth=1.5;
            ax=H.CurrentAxes;ax.FontName='Calibri';ax.Box='off';
            ax.FontSize=NSPF.atlSize;ax.FontWeight='bold';
            %ax.XTick=[1,2];ax.XTickLabel=xTickLabels;
            ax.YLabel.String = 'Pearson Correlation';
            ax.YLabel.FontSize = NSPF.alSize;ax.YLabel.FontWeight = 'bold';
            ax.XLabel.String = 'Neuron #';
            ax.YLabel.FontSize = NSPF.alSize;ax.YLabel.FontWeight = 'bold';
            %ax.XLim=[0,0.3];ax.YLim=[-inf,inf];
            ax.LineWidth=2;
        end
        function PlotDFF_RJONOFF(dFF_raw,RJONOFF,vel_t,vel,RJVel)
            spacing=std(dFF_raw(:))*7;
            [Ntrace,nT]=size(dFF_raw);
            Height=(Ntrace+2)*spacing;
            sum_dFF=sum(dFF_raw,1);
            figure('Name', 'Raw DFF w/ RJ ON/OFF overlay');
            hold on;
            % draw RJON/OFF color overlay
            for i=1:numel(RJONOFF)
                start=1000*(i-1)+1;
                stop=1000*i;
                p1=[start,-spacing];p2=[stop,-spacing];p3=[stop,Height];p4=[start,Height];
                XCoord=[p1(1),p2(1),p3(1),p4(1)];
                YCoord=[p1(2),p2(2),p3(2),p4(2)];
                if (RJONOFF(i))
                    patch(XCoord,YCoord,[0.4660 0.6740 0.1880],'EdgeColor','none')
                else
                    patch(XCoord,YCoord,[0.9290 0.6940 0.1250],'EdgeColor','none')
                end
            end
            %draw dF/F
            for i=1:size(dFF_raw,1)
                plot(dFF_raw(i,:)+(size(dFF_raw,1)-i)*spacing,'k','LineWidth',2);
            end
            % draw sum DFF
            plot((sum_dFF/max(sum_dFF)*spacing-spacing)*3,'r');
            % draw video velocity
            plot(vel_t,(vel/50-2)*spacing*3,'b');
            % draw RJ Angle_vel
            plot((RJVel/max(RJVel)-3)*spacing*3,'b')
            % turn off axes
            set(gca,'visible','off');
        end
        function filtered=FilterData(raw,filter_vec)
            [N,~]=size(raw);
            filtered=raw;
            filter_mat=repmat(filter_vec,[N,1]);
            filtered(filter_mat==0)=[];
            filtered=reshape(filtered,[N,numel(filtered)/N]);
        end
        function PlotDFF_Selected(dFF_raw,ind,behavior,figTitle)
            spacing=std(dFF_raw(:))*7;
            figure('Name', figTitle);
            hold on;
            %draw dF/F
            for i=1:numel(ind)
                ii=ind(i);
                plot(dFF_raw(ii,:)+(numel(ind)-i+1)*spacing,'k','LineWidth',2);
            end
            % draw sum DFF
            plot(behavior/std(behavior)/5*spacing,'b');
%             if i==1
%                 spikeTrace=dFF_raw(ii,:);
%                 peakTh=0.2*max(spikeTrace);
%                 [PKs,LOCs]=findpeaks(spikeTrace,'MinPeakHeight',peakTh,'MinPeakWidth',3);
%                 plot(LOCs,PKs+spacing,'xb');
%             end
            title(sprintf('dFF spacing=%.2f   behavior spacing=%.2f',spacing,5*std(behavior)));
            axis off;
            %set(gca,'visible','off');
        end
        function PlotDFF_Sessioned(dFF_raw,behavior,overlay,sessionLen,title,color)
            if nargin==5
                C=colororder;
                C=repmat(C,[10,1]);
            else
                C=color;
            end
            [N,nT]=size(dFF_raw);
            nS=round(nT/sessionLen);
            spacing=std(dFF_raw(:))*7;
            figure('Name', title);
            hold on;
            %draw overlay
            
            idx_overlay=find(overlay==1);
            X=repmat(idx_overlay,[2,1]);
            Y=repmat([0;(N+3)*spacing],[1,numel(idx_overlay)]);
            line(X,Y,'Color',[0.7,0.7,0.7]);
            %draw dF/F
            for i=1:N
                plot(dFF_raw(i,:)+(N-i+2)*spacing,'k','LineWidth',2);
            end
            % draw sessioned behavior
            
            scale_beh=1/std(behavior)/5*spacing*2;
            for i=1:nS
                idx=((i-1)*sessionLen+1):i*sessionLen;
                plot(idx,behavior(idx)*scale_beh,'Color',C(i,:));
                if i~=nS
                    line([max(idx),max(idx)],[0,(N+3)*spacing],...
                        'Color','blue','LineWidth',2,'LineStyle','--');
                end
            end
            set(gca,'visible','off'); 
        end
        function PlotDFF_Colored(dFF_raw,scale,colorList,title)
            figure('Name', title);
            hold on;
            [N,~]=size(dFF_raw);
            for i=1:N
                plot(dFF_raw(i,:)+(N-i+1)*scale,'LineWidth',2,'Color',colorList(i,:));
            end
            set(gca,'visible','off'); 
        end
        function I=PlotMap_Highlight(Coor,size,upsample,ind,theme)
            if nargin==4
                theme=1;
            end
            if theme==1
                I=zeros(size*upsample,size*upsample);
            else
                if theme==2
                    I=ones(size*upsample,size*upsample);
                    I=insertShape(I,'FilledCircle',...
                        [size*upsample/2,size*upsample/2,size*upsample/2],...
                        'Color',[0.7,0.7,0.7],'LineWidth',1);
                end
            end
            for i=1:numel(Coor)
                coord=Coor{i}*upsample;
                %centroid=mean(coord,2);
                coord=coord(:);
                if sum(i==ind)>0
                    color='red';
                else
                    if theme==1
                        color='white';
                    else
                        if theme==2
                            color='black';
                        end
                    end
                end
                I=insertShape(I,'Polygon',coord','Color',color,'LineWidth',2);
                %I=insertShape(I,'FilledPolygon',coord','Color','white','Opacity',0.2);
                %I=insertText(I,centroid',sprintf('%d',i),'TextColor','white',...
                %    'AnchorPoint','Center','BoxOpacity',0,'FontSize',24);
            end
            figure('Name','Neuron Map(Highlighted)');
            image(I);
            axis image;
        end
        function [avg,n]=GetSpikeTriggerAvg(spikeTrace,sessionLen,behavior,hlfW,nAvg,scale)
            if nargin==4
                nAvg=1;
                scale=-1;
            else
                if nargin==5
                    scale=-1;
                end
            end
            nS=numel(spikeTrace)/sessionLen;
            peakTh=0.2*max(spikeTrace);
            traceIdx=1:numel(spikeTrace);
            if scale~=-1
                behavior=behavior*scale;
            end
            behaviorIdx=(1:numel(behavior))/numel(behavior)*numel(traceIdx);
            multiplier=numel(behaviorIdx)/numel(traceIdx);
            hlfW_beh=ceil(hlfW*multiplier);
            epochs=[];
            for i=1:nS
                sessionIdx=(i-1)*sessionLen+1:i*sessionLen;
                spikeTraceSession=spikeTrace(sessionIdx);
                if sum(spikeTraceSession>=peakTh)<1
                    continue;
                end
                behSessionIdx=round(sessionIdx*multiplier);
                behSession=behavior(behSessionIdx);
                behSessionPadded=[zeros(1,hlfW_beh),behSession,zeros(1,hlfW_beh+1)];                
                [~,locsSession]=findpeaks(spikeTraceSession,'MinPeakHeight',peakTh,'MinPeakWidth',3);
                for j=1:numel(locsSession)
                    begin=locsSession(j)-hlfW_beh+hlfW_beh;
                    finish=begin+2*hlfW_beh;
                    epoch=behSessionPadded(begin:finish);
                    epoch=movmean(epoch,nAvg);
                    epochs=[epochs;epoch];
                end
            end
            avg=mean(epochs,1);
            sigma=std(epochs,1);
            [n,~]=size(epochs);
            epoch_t=(-hlfW_beh:hlfW_beh)/multiplier/3;
            hold on;
            %plot(epoch_t,epochs,'Color',[0.9 0.9 0.9]);
            pAvg=plot(epoch_t,avg,'b','LineWidth',3);
            pStd=plot(epoch_t,sigma,'Color',[0.7,0.7,0.7]);
            %plot(epoch_t,movmean(avg,8),'b','LineWidth',2);
            ymin=min(avg)*0.9; ymax=max(avg)*1.1;
            axis([min(epoch_t),max(epoch_t),ymin,ymax])
            line([0;0],[ymin;ymax],'LineWidth',1,'Color','black','LineStyle','--');
            if scale~=-1
                ax=gca;ax.LineWidth=2;ax.FontName='Calibri';
                ax.FontSize=NSPF.atlSize;ax.FontWeight='bold';
                ax.YLabel.String = 'Locomition Speed (mm/s)';
                ax.YLabel.FontSize =NSPF.alSize;ax.YLabel.FontWeight = 'bold';
                ax.YLim=[ymin,ymax];
                ax.XLabel.String = 'Time (s)';
                ax.YLabel.FontSize =NSPF.alSize;ax.YLabel.FontWeight = 'bold';
                l=legend([pAvg,pStd],'Mean','STD');
                l.Box='off';
            end
        end
        function [TrajIm,distance]=PlotTrajectories(trajectories,flags,size,figTitle,color,hFigure)
            if nargin==4
                color='b';
                hFigure=figure('Name',figTitle);
            end
            hold on;
            th=NSPF.speedTh/400*size;
            distance=0;
            for i=1:numel(trajectories)
                if (flags(i))
                    trajData=trajectories{i};
                    traj_interp=NSPF.TrajInterp(trajData(1),0.95,20);
                    x=traj_interp.x;
                    y=-traj_interp.y;
                    plot(x+1i*y,'Color',color);
                    dx=[0,diff(x)];
                    dy=[0,diff(y)];
                    del=sqrt(dx.^2+dy.^2);
                    del(del<th)=0;
                    distance=distance+sum(del)/size;
                end
            end
            title(sprintf('Raw distance traversed: %.3f',distance))
            axis equal;
            axis([0,size,-size,0]);
            axis off;
            F = getframe(hFigure);
            [TrajIm, ~] = frame2im(F);
        end
        function [TrajIm,distance]=PlotTrajectory(trajectory,range,size,figTitle,color,hFigure)
            if nargin==4
                color='b';
                hFigure=figure('Name',figTitle);
            end
            hold on;
            th=NSPF.speedTh/400*size;



            trajData=trajectory;
            traj_interp=NSPF.TrajInterp(trajData(1),0.95,20);
            x=traj_interp.x(range);
            y=-traj_interp.y(range);
            plot(x+1i*y,'Color',color);
            dx=[0,diff(x)];
            dy=[0,diff(y)];
            del=sqrt(dx.^2+dy.^2);
            del(del<th)=0;

            axis equal;
            axis([0,size,-size,0]);
            axis off;
            F = getframe(hFigure);
            [TrajIm, ~] = frame2im(F);
        end
        function stats=BehaviorStats(behavior,nVec,RJONOFF,th)
            stats.totDist=zeros(1,sum(RJONOFF));
            stats.avgSpeed=zeros(1,sum(RJONOFF));
            startIdx=1;
            ii=0;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                if RJONOFF(i)
                    ii=ii+1;
                    epoch=behavior(startIdx:endIdx);
                    stats.totDist(ii)=sum(epoch(epoch>th));
                    stats.avgSpeed(ii)=sum(epoch(epoch>th))/sum(epoch>th);
                end
                startIdx=startIdx+nVec(i);
            end
        end
        function stats=AngleStats(angle,nVec,RJONOFF)
            stats.totAngle=zeros(1,sum(RJONOFF));
            startIdx=1;
            ii=0;
            for i=1:numel(nVec)
                endIdx=startIdx+nVec(i)-1;
                if RJONOFF(i)
                    ii=ii+1;
                    epoch=unwrap(angle(startIdx:endIdx));
                    %stats.totAngle(ii)=abs(epoch(end)-epoch(1));
                    diff_epoch=abs(diff(epoch));
                    diff_epoch(diff_epoch<0.03)=0;
                    stats.totAngle(ii)=sum(diff_epoch);
                end
                startIdx=startIdx+nVec(i);
            end
        end     
        function Rates=PlotFluoRates(dFF_raw,flag,FPS,xTickLabels)
            Rates.Frates1=sum(dFF_raw(:,flag),2)/sum(flag)*FPS;
            Rates.Frates2=sum(dFF_raw(:,flag==0),2)/sum(flag==0)*FPS;
            meanRates=[mean(Rates.Frates1),mean(Rates.Frates2)];
            Rates.meanRates=meanRates;
            stdRates=[std(Rates.Frates1),std(Rates.Frates2)];
            Rates.stdRates=stdRates;
            H=figure('Name','Fireing Rates');
            hold on;
            p=plot([Rates.Frates1';Rates.Frates2'],'Color',[0.7,0.7,0.7],'LineWidth',1);
            pMean=plot(meanRates,'k','LineWidth',3);
            ax=H.CurrentAxes;ax.FontName='Calibri';
            ax.FontSize=NSPF.atlSize;
            ax.YAxis.FontSize=NSPF.atlSize;
            ax.XAxis.FontSize=NSPF.alSize;ax.FontWeight='bold';
            ax.XTick=[1,2];ax.XTickLabel=xTickLabels;
            ax.YLabel.String = 'Fluorescence rates (AU)';
            ax.YLabel.FontSize =NSPF.alSize;ax.YLabel.FontWeight = 'bold';
            ax.XLim=[0.75,2.25];ax.YLim=[-inf,inf];
            ax.LineWidth=2;
            l=legend([p(1),pMean],'Individual','Mean');
            l.Box='off';
            ax.PlotBoxAspectRatio=[1,1.5,1];
        end
        function PlotBestCorrelatedSession(dFFTrace,behavior,sessionLen)
            nS=numel(dFFTrace)/sessionLen;
            corr=zeros(1,nS);
            p=zeros(1,nS);
            for i=1:nS
                idx=(i-1)*sessionLen+1:i*sessionLen;
                dFFSession=dFFTrace(idx);
                behSession=behavior(idx);
                [mat,pmat]=corrcoef(dFFSession,behSession);
                corr(i)=mat(1,2);
                p(i)=pmat(1,2);
            end
            corr(p>0.05)=-inf;
            ii=find(corr==max(corr));
            idx=(ii-1)*sessionLen+1:ii*sessionLen;
            spacing=std(dFFTrace(idx))*7;
            figure('Name', 'Highest correlated dFF session');
            hold on;
            plot(dFFTrace(idx)+spacing,'k','LineWidth',2);
            plot(behavior(idx)/std(behavior(idx))/5*spacing,'b');
            set(gca,'visible','off');
        end
        function H=MyBar(data,xTickLabels,xLabel,yLabel,legendText,axLim,AR,title)
            H=figure('Name',title);
            hold on;
            b=bar(data);
            for i=1:numel(b)
                b(i).LineWidth=1.5;
            end
            ax=H.CurrentAxes;ax.FontName='Calibri';
            ax.FontSize=NSPF.atlSize;
            ax.YAxis.FontSize=NSPF.atlSize;
            ax.XAxis.FontSize=NSPF.alSize;ax.FontWeight='bold';
            ax.XTickLabel=xTickLabels;
            ax.XLabel.String = xLabel;
            ax.XLabel.FontSize =NSPF.alSize;ax.XLabel.FontWeight = 'bold';
            ax.YLabel.String = yLabel;
            ax.YLabel.FontSize =NSPF.alSize;ax.YLabel.FontWeight = 'bold';
            ax.XLim=axLim(1,:);ax.YLim=axLim(2,:);
            ax.LineWidth=2;
            l=legend(legendText);
            l.Box='off';
            ax.PlotBoxAspectRatio=[AR,1];
        end
        function pval=PValCalc(X,Y)
            %pval=0.5+0.5*erf(-mean(X-Y)/sqrt(2)/sqrt(std(X)^2/numel(X)+std(Y)^2/numel(Y)));
            pval=0.5+0.5*erf(-mean(X-Y)/sqrt(2)/sqrt(std(X)^2/numel(1)+std(Y)^2/numel(1)));
        end
    end
end