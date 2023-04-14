
% this script needs the original pclamp (Molecular devices) .abf files to
% run because collecting all the recordings into a single matlab file
% results in a dataset too large to save to GitHub.  If you would like the
% original data set please contact the John H. Byrne or Curtis L. Neveu.



folder = 'D:\Data\2019\';
[info1,textr] = xlsread([folder 'data_info.xlsx'],'gapfree');


continuous = false; % plots the entire window
plateau_w = [1,2];% sets window for detecting plateau (s)
Istep = true; % measures current as change from baseline
bl_w = [-0.2,0]; % window for measuring baseline current (s)
plotindIO = false; % plot individual IO traces 
plt_res = false; % plot input resistance trace
plt_res_all = false; % plot all resistance traces
irW = [-0.1,0]; % time to measure input resistance
pthr = -70; % threshold for plateau

correlation = false; % plot the correlation of voltage between electrodes
corrW = 5;
c_step = 2000;

ST_trace = false; % plot overlapped traces during step

PXpum10 = 0.976; %pixels per um for 10X


info1(:,14:15) = info1(:,10:11)*PXpum10; 
info1(:,16) = (info1(:,15)./4).^2;% calculate area
info1(:,17) = (info1(:,15) - info1(:,14).*2)./2; % average axon width

if plotindIO
    close(findobj(0,'Name', 'IOi'))
    figIOi = figure('Name','IOi','NumberTitle','off');
end

cnt = true;
sp_cnt = cell(size(info1,1),3);
sp_i = sp_cnt;
meanI = 0.1:0.1:4;
if ~Istep; meanI = -4:0.1:2;end
meansp = nan([size(meanI,2),size(sp_i,1),3]);
plateaum = false(size(meansp));
plateauI = nan(size(info1,1),3);
plateau_dur = nan(size(info1,1),3);
IR = plateauI;
ihold = nan(size(info1,1),3);
if plt_res_all
    close(findobj(0,'Name', 'IRtracesall'))
    figure('Name','IRtracesall','NumberTitle','off');
end
for f=1:size(info1,1)
    [data,si,h] = abfload([folder,num2str(info1(f,1)),'.abf']);
    if f==2; data(:,[1 3]) = data(:,[1 3])*10;end  % for some reason 2nd data is off by factor of 10
    if cnt 
        info1(:,2:9) = round(info1(:,2:9)*1000/si);
        X = 0:si/1e6:size(data,1)*si/1e6;
    end
    
    for w=1:3
        if ~isnan(info1(f,(w-1)*2+8))
            ch = 2;
            if f>32; ch = 3; end
            rdata = data(info1(f,(w-1)*2+8):info1(f,(w-1)*2+9),ch:end);
            ron = find(diff(rdata(:,2))<-0.4,1);
            roff = find(diff(rdata(:,2))>0.4,1);
            if f==1
                rdata(ron:roff,1) = rdata(ron:roff,1)+9;% the current was accidently injected into ME1 for this exp (corrects electrode voltage drop)
            end
            V = mean(rdata(roff+irW(1)*5000:roff+irW(2)*5000,1));
            I = mean(rdata(roff+irW(1)*5000:roff+irW(2)*5000,2));
            IR(f,w) = V/I;

            if plt_res
                close(findobj(0,'Name', 'IRtraces'))
                figure('Name','IRtraces','NumberTitle','off');

                axr(1) = subplot(2,1,1);
                 plot(X(1:size(rdata,1)),rdata(:,1));hold on
                 plot(X([roff+irW(1)*5000, roff+irW(2)*5000]),[IR(f,w),IR(f,w)],'r','LineWidth',2)
                axr(2) = subplot(2,1,2);
                 plot(X(1:size(rdata,1)),rdata(:,2));hold on  
                linkaxes(axr(:),'x')
            end
        

            excp = data(info1(f,(w-1)*2+2):info1(f,(w-1)*2+3),ch:end);

%             if w==1% measuring the holding current
            ihold(f,w) = mean(excp(1:15000,2));
%             end

            spike = find(excp(:,1)>0);
            spike(find(diff(spike)<5)+1) = [];% remove duplicate detection

            step = 5;
            on = find(diff(excp(1:step:end,2))>0.09);
            on(find(diff(on)<5)+1) = [];% remove duplicate detection
            off = find(diff(excp(1:step:end,2))<-0.09);
            off(find(diff(off)<5)+1) = [];% remove duplicate detection
            on = on*step;
            off = off*step;
            misson = arrayfun(@(x) ~any(abs(on-x+5000)<100),off);
            on = sort([on;off(misson)-5000]);% add to on if off detected stim but on did not
            missoff = arrayfun(@(x) ~any(abs(off-x-5000)<100),on);
            off = sort([off;on(missoff)+5000]);

            logg = off-on<5000;
            on(logg) = [];
            off(logg) = [];

            sp_cnt{f,w} = zeros(size(on));
            sp_i{f,w} = zeros(size(on));
            plateau = false(size(on));
            spike_m_log = false(size(spike));
            
            for s=1:length(on)
                sp_cnt{f,w}(s) = sum(spike>on(s) & spike<off(s));
                sp_i{f,w}(s) = mean(excp(on(s):off(s),2));
                if Istep
                    sp_i{f,w}(s) = sp_i{f,w}(s) - mean( excp(on(s)+bl_w(1)*1e6/si : on(s)+bl_w(2)*1e6/si, 2) );
                end
                sp_i{f,w}(s) = round(sp_i{f,w}(s)*10)/10;
                plateau(s) = mean(excp(off(s)+plateau_w(1)*1e6/si:off(s)+plateau_w(2)*1e6/si,1))>pthr;
                if ~plateau(s)
                    spike_m_log(spike>on(s) & spike<off(s)) = true;
                end
            end
            if f==33 && w==1
                plateau(:) = false;% a spurious depolarization, no actual plateaus present.
            end
            plateaum(abs(meanI - sp_i{f,w}(plateau))<0.001,f,w) = true;
            if any(plateau)
                plateauI(f,w) = meanI(abs(meanI - sp_i{f,w}(find(plateau,1)))<0.001);
                if find(plateau)<length(plateau)
                    idx = find(excp(on(plateau):on(find(plateau)+1),1)>pthr)+on(plateau);
                else
                    idx = find(excp(on(plateau):end,1)>pthr)+on(plateau);
                end
                idx(idx<on(plateau)+5000) = [];
                plateau_dur(f,w) = length(idx)/5000;
            end

            sp_cnt{f,w}(plateau) = [];
            sp_i{f,w}(plateau) = [];
            on(plateau) = [];
            off(plateau) = [];
            spike_m = spike(spike_m_log);

            stimW = [-0.1,1.1]; % window for the stimuli
            st_data = zeros(ceil(diff(stimW)*1e6/si), length(on), 2);
            for s=1:length(on) % combine all the stimuli to one graph
                st_data(:,s,:) = excp(on(s)+stimW(1)*1e6/si:on(s)+stimW(2)*1e6/si,1:2);
            end

            if ST_trace
                close(findobj(0,'Name', 'IOtrace'))
                figure('Name','IOtrace','NumberTitle','off');
                sax(1) = subplot(2,1,1);
                 Xs = repmat(linspace(stimW(1),stimW(2),size(st_data,1))',1,length(on));
                 plot(Xs,st_data(:,:,1))
                sax(2) = subplot(2,1,2);
                 plot(Xs,st_data(:,:,2))
                linkaxes(sax(:),'x')
            end

            for s=1:size(sp_i{f,w},1)
                logicc = abs(sp_i{f,w}-sp_i{f,w}(s))<0.001;
                meansp(abs(meanI - sp_i{f,w}(s))<0.001,f,w) = mean(sp_cnt{f,w}(logicc));
            end

            txt = regexprep(split(mat2str(1:size(sp_i{f,w},1)),' '),'\D','');
            
            if plt_res_all && w==1
                figure(findobj(0,'Name', 'IRtracesall'))
                axpp = subplot(5, ceil(size(info1,1)/5), f);
                cclr = 'rb';
                 plot(X(ron-4000:roff+5000),rdata(ron-4000:roff+5000,1),cclr(any(plateau)+1));hold on               
                axpp.YLim = [-140,-75];
            end

            if continuous && any(plateau)% && f>32%&& f>3 && f<12%&& 
                np = 2;
                if correlation; np = 3; end

%                 close(findobj(0,'Name', 'traces'))
%                 figure('Name','traces','NumberTitle','off');
                figure('Name',[num2str(f) '  ' num2str(info1(f,1)) '  '  num2str(w)],'NumberTitle','off');

                ax(1) = subplot(np,1,1);
                 plot(X(1:size(excp,1)),excp(:,1));hold on
                 scatter(X(spike_m),excp(spike_m,1),'d');
                 plot(idx/5000,repelem(pthr,length(idx)));hold on;
                ax(2) = subplot(np,1,2);
                 plot(X(1:size(excp,1)),excp(:,2));hold on
    %              fth = diff(excp(1:step:end,2));  
    %              plot(X(1:length(fth))*step,fth);hold on
                 M = mean(excp(:,2));
                 scatter(X(on),ones(length(on),1)*M ); hold on
                 scatter(X(off),ones(length(off),1)*M);hold on
                 if any(plateau)
                    plot( [X(on);X(off)], ones(2,length(on))*M, 'LineWidth',3); hold on
                 end
                 text( mean([X(on);X(off)]), ones(1,length(on))*M, txt,...
                     'HorizontalAlignment','center','VerticalAlignment','bottom')      
                
                
                if correlation
                    MEcorr = zeros(floor(length(excp)/c_step),1);% measure the correlation between electrodes
                    y = downsample(excp,c_step);
                    for i=1:length(MEcorr)-corrW
                        dta = y(i:i+corrW,[1 3]);
                        MEcorr(i) = corr(dta(:,1),dta(:,2));
                    end   

                    ax(3) = subplot(np,1,3);
                     plot(X(1:c_step:length(MEcorr)*c_step),MEcorr)
                end
                linkaxes(ax(:),'x')
            end
        end
    end    
    
    if plotindIO
        figure(figIOi)
        plot(sp_i{f,1},sp_cnt{f,1},'k-');hold on
        plot(sp_i{f,2},sp_cnt{f,2},'r-');hold on
        %text(sp_i,sp_cnt+1,txt,'HorizontalAlignment','center','VerticalAlignment','bottom')
    end
    cnt = false;
end

%%
% close(findobj(0,'Name', 'IO'))
% figIO = figure('Name','IO','NumberTitle','off');
% 
% sterr = squeeze(nanstd(meansp,0,2)/sqrt(size(meansp,2)) );
% msp = squeeze(nanmean(meansp,2));
% 
% nmeanI = meanI;
% nmeanI(any(isnan(msp),2))=[];
% sterr(any(isnan(sterr),2),:)=[];
% msp(any(isnan(msp),2),:)=[];
% 
% color = 'kr';
% for p=1:2
%     x_plot = [nmeanI fliplr(nmeanI)];                                         % Transposed To Row Vectors
%     y_plot = [msp(:,p)'-sterr(:,p)' fliplr(msp(:,p)'+sterr(:,p)')];       % Transposed To Row Vectors
%     plt(p) = plot(nmeanI,msp(:,p),color(p),'LineWidth',3);hold on %#ok<SAGROW>
%     fill(x_plot, y_plot, color(p), 'FaceAlpha',0.5,'EdgeColor',color(p));hold on
% end
% legend(plt,'pretest','Ni')
% axx = gca;
% axx.YLabel.String = 'Number of spikes';
% axx.XLabel.String = 'Current step (nA)';
% 
% 
% close(findobj(0,'Name', 'Axon'))
% figIO = figure('Name','Axon','NumberTitle','off');
% 
% 
% subplot(2,3,1)
% scatter(info1(:,10),ihold)
% axss(1) = gca;
% axss(1).XLabel.String = 'Axon length (\mum)';
% axss(1).YLabel.String = 'Holding current (nA)';
% 
% 
% subplot(2,3,2)
% scatter(info1(:,10),plateauI(:,2))
% axss(3) = gca;
% axss(3).XLabel.String = 'Axon length (\mum)';
% axss(3).YLabel.String = 'Plateau threshold (nA)';
% 
% subplot(2,3,3)
% scatter(info1(:,10),IR(:,1))
% axss(4) = gca;
% axss(4).XLabel.String = 'Axon length (\mum)';
% axss(4).YLabel.String = 'Input Resistance (M\Omega)';
% 
% 
% pp = any(plateaum(:,:,1));
% group = repelem("",length(ihold));
% group1 = group;
% group(pp) = 'Plateau';
% group(~pp) = 'None';
% 
% disp(' ')
% disp('Axon length of plateau vs none')
% %ax = barplotit(info(4:end,10),group(4:end),subplot(2,3,4));
% ax = barplotit(info1(:,10),group,subplot(2,3,4));
% ax.YLabel.String = 'Axon length (\mum)';
% 
% disp('Input resistance of plateau vs none')
% ax = barplotit(IR(:,1),group,subplot(2,3,5));
% ax.YLabel.String = 'Input Resistance (M\Omega)';
% 
% disp('Holding current of plateau vs none')
% ax = barplotit(ihold(:,1),group,subplot(2,3,6));
% ax.YLabel.String = 'Holding Current (nA)';

%%  Miami vs South Coast

% close(findobj(0,'Name', 'Source'))
% figIO = figure('Name','Source','NumberTitle','off');
% 
% miami = cellfun(@(s) contains(s,'M'),textr(2:end,1))';
% 
% subplot(2,4,1)
% mp = sum(miami & pp);    ma = sum(miami);
% scp = sum(~miami & pp);  sca = sum(~miami);
% source = [repmat('m',ma,1); repmat('s',sca,1)];
% plat = [ones(mp,1); zeros(ma-mp,1);ones(scp,1);zeros(sca-scp,1)];
% [tbl,chi2stat,pval] = crosstab(source,plat);
% fprintf('\nUncut vs Cut:  Counts: Miami (%i;%i) SC (%i;%i) X = %.1f , p = %.3f\n',mp,ma,scp,sca,chi2stat,pval);
% bar([mp/ma scp/sca]*100)
% ax = gca;
% ax.XTickLabel = {'Miami','SC'};
% ax.YLabel.String = 'Percent Plateau';
% ytickformat(ax, 'percentage')
% 
% jit = (rand(size(ihold'))-0.5)*.2;
% 
% subplot(2,4,2)
% scatter(~miami+1+jit,ihold,15)
% men =  prctile(ihold(miami),[25 50 75]);
% rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
% men2 =  prctile(ihold(~miami),[25 50 75]);
% rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
% line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
% ax = gca;
% ax.XLim = [0,3];
% ax.XTick = [1,2];
% ax.XTickLabels = {'Miami','SC'};
% ax.YLabel.String = 'Holding Current (nA)';
% 
% subplot(2,4,3)
% scatter(~miami+1+jit,IR(:,1),15)
% men =  prctile(IR(miami,1),[25 50 75]);
% rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
% men2 =  prctile(IR(~miami,1),[25 50 75]);
% rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
% line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
% ax = gca;
% ax.XLim = [0,3];
% ax.XTick = [1,2];
% ax.XTickLabels = {'Miami','SC'};
% ax.YLabel.String = 'Input Resistance (M\Omega)';
% 
% 
% 
% 
% uncut = cellfun(@(s) contains(s,'A'),textr(2:end,1))';
% 
% subplot(2,4,5)
% mp = sum(uncut & pp);    ma = sum(uncut);
% scp = sum(~uncut & pp);  sca = sum(~uncut);
% source = [repmat('m',ma,1); repmat('s',sca,1)];
% plat = [ones(mp,1); zeros(ma-mp,1);ones(scp,1);zeros(sca-scp,1)];
% [tbl,chi2stat,pval] = crosstab(source,plat);
% fprintf('Uncut vs Cut:  Counts: Uncut (%i;%i) Cut (%i;%i) X = %.1f , p = %.3f\n',mp,ma,scp,sca,chi2stat,pval);
% bar([mp/ma scp/sca]*100)
% ax = gca;
% ax.XTickLabel = {'Uncut','Cut'};
% ax.YLabel.String = 'Percent Plateau';
% ytickformat(ax, 'percentage')
% 
% subplot(2,4,6)
% jit = (rand(size(uncut))-0.5)/10;
% scatter(~uncut+1+jit,ihold,15)
% men =  prctile(ihold(uncut),[25 50 75]);
% rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
% men2 =  prctile(ihold(~uncut),[25 50 75]);
% rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
% line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
% ax = gca;
% ax.XLim = [0,3];
% ax.XTick = [1,2];
% ax.XTickLabels = {'Uncut','Cut'};
% ax.YLabel.String = 'Holding Current (nA)';
% 
% subplot(2,4,7)
% scatter(~uncut+1+jit,IR(:,1),15)
% men =  prctile(IR(uncut,1),[25 50 75]);
% rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
% men2 =  prctile(IR(~uncut,1),[25 50 75]);
% rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
% line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
% ax = gca;
% ax.XLim = [0,3];
% ax.XTick = [1,2];
% ax.XTickLabels = {'Uncut','Cut'};
% ax.YLabel.String = 'Input Resistance (M\Omega)';
% 
% 
% subplot(2,4,8)
% scatter(~uncut+1+jit,info1(:,10),15)
% men =  prctile(info1(uncut,10),[25 50 75]);
% rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
% men2 =  prctile(info1(~uncut,10),[25 50 75]);
% rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
% line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
% ax = gca;
% ax.XLim = [0,3];
% ax.XTick = [1,2];
% ax.XTickLabels = {'Uncut','Cut'};
% ax.YLabel.String = 'Axon length (\mum)';
% % yaxis = {'Axon length (\mum)','Axon circumferance (\mum)','Axon area','Axon width (\mum)'};
% % for s=1:4
% %     subplot(1,4,s)
% %     scatter(ihold,info(:,9+s))
% %     axss = gca;
% %     axss.YLabel.String = yaxis{s};
% %     axss.XLabel.String = 'Holding potential (nA)';
% % end

%%

close(findobj(0,'Name', 'Plateau Duration'))
figIO = figure('Name','Plateau Duration','NumberTitle','off');

plateau_dur(isnan(plateau_dur(:,1)),1) = 0;
Y = plateau_dur(4:11,1:2);


subplot(1,2,1)
men =  prctile(Y(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(Y(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(Y,1)),Y')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+}'};
ax.YLabel.String = 'Duration (s)';

dY = diff(Y,[],2);
[hn, pnorm, W] = adtest(dY);
[h,p,ci,stats] = ttest(dY);
[coef, pval] = corr(Y);
disp(['pvalue = ' num2str(p)])
disp(stats)

plateau_dur(isnan(plateau_dur(:,1)),1) = 0;
Y = plateau_dur(33:end,[1 3]);


subplot(1,2,2)
men =  prctile(Y(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(Y(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(Y,1)),Y')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+} + Cd^{2+}'};
ax.YLabel.String = 'Duration (s)';

dY = diff(Y,[],2);
[hn, pnorm, W] = adtest(dY);
[h,p,ci,stats] = ttest(dY);
[coef, pval] = corr(Y);
%% spike threshold

threshold_sp = sp_cnt(:,1:2);
threshold_i = sp_i(:,1:2);

rmv = any(cellfun(@isempty,threshold_sp),2);
threshold_sp(rmv,:) = [];
threshold_i(rmv,:) = [];

thridx = cellfun(@(x) find(x,1),threshold_sp);
thr = cellfun(@(x,y) x(y),threshold_i,num2cell(thridx));

close(findobj(0,'Name', 'Threshold'))
figIO = figure('Name','Threshold','NumberTitle','off');

subplot(1,2,1)
men =  prctile(thr(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(thr(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(thr,1)),thr')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+}'};
ax.YLabel.String = 'Spike threshold (nA)';

dYthr{1} = diff(thr,[],2);
[hnthr, pnormthr, Wthr] = adtest(dYthr{1});
[hthr,pthr,cithr,statsthr] = ttest(dYthr{1});

[coef, pval] = corr(Y);


threshold_sp = sp_cnt(:,[1 3]);
threshold_i = sp_i(:,[1 3]);

rmv = any(cellfun(@isempty,threshold_sp),2);
threshold_sp(rmv,:) = [];
threshold_i(rmv,:) = [];

thridx = cellfun(@(x) find(x,1),threshold_sp);
thr = cellfun(@(x,y) x(y),threshold_i,num2cell(thridx));

subplot(1,2,2)
men =  prctile(thr(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(thr(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(thr,1)),thr')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+}'};
ax.YLabel.String = 'Spike threshold (nA)';

dYthr{2} = diff(thr,[],2);
[hnthr, pnormthr, Wthr] = adtest(dYthr{2});
[hthr,pthr,cithr,statsthr] = ttest(dYthr{2});
%% holding

close(findobj(0,'Name', 'holding'))
figIO = figure('Name','holding','NumberTitle','off');

nihold = ihold(~any(isnan(ihold(:,1:2)),2),1:2);

subplot(1,2,1)
men =  prctile(nihold(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(nihold(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(nihold,1)),nihold')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+}'};
ax.YLabel.String = 'holding current (nA)';

dYhold{1} = diff(nihold,[],2);
[hnhold, pnormhold, Whold] = adtest(dYhold{1});
[hhold,phold,cihold,statshold] = ttest(dYhold{1});

% [hcoef, hpval] = corr(dYthr,dYhold);


nihold = ihold(~any(isnan(ihold(:,[1 3])),2),[1 3]);

subplot(1,2,2)
men =  prctile(nihold(:,1),[25 50 75]);
rectangle('Position',[0.8, men(1),0.4,men(3)-men(1)]);hold on
men2 =  prctile(nihold(:,2),[25 50 75]);
rectangle('Position',[1.8, men2(1),0.4,men2(3)-men2(1)]);hold on
line([0.75 1.75;1.25 2.25],[men(2) men2(2);men(2) men2(2)],'Color','k','LineWidth',2)
plot(repmat((1:2)',1,size(nihold,1)),nihold')
ax = gca;
ax.XLim = [0,3];
ax.XTick = [1,2];
ax.XTickLabels = {'Control','Ni^{2+}'};
ax.YLabel.String = 'holding current (nA)';

dYhold{2} = diff(nihold,[],2);
[hnhold, pnormhold, Whold] = adtest(dYhold{2});
[hhold,phold,cihold,statshold] = ttest(dYhold{2});

[hcoef, hpval] = corr(dYthr{2},dYhold{2});