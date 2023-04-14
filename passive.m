
% this script needs the original pclamp (Molecular devices) .abf files to
% run because collecting all the recordings into a single matlab file
% results in a dataset too large to save to GitHub.  If you would like the
% original data set please contact the John H. Byrne or Curtis L. Neveu.

folder = 'D:\Data\2019\';

file = ["Atype","Delayed_K","KCa","NaP","Ltype","Rtype"];
neuron = ["B51","B64","B8"];

opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
pass = zeros(0,6);
group = cell(0,1);
cfun = @(p,x) p(1).*exp(-x./p(2)) + p(3);
cp0 = [20,0.2,-80];
info = cell(3);
cnt = 1;
nn = [3,2,3,3,3,3];
datac = cell(3);
for f=1:length(file)
    fprintf(['\n\n\n>>>>>' file{f} '<<<<<<<\n\n\n'])
    for n=1:nn(f)
        [~,~,info{n,f}] = xlsread([folder file{f} '.xlsx'],neuron{n});
        info{n,f} = info{n,f}(3:end,2:end);
        info{n,f} = cellfun(@num2str,info{n,f},'UniformOutput',false);
        
        if strcmp(file{f},'NaP')
            datac{n,f} = readpass(folder,info{n,f}(:,5:7));
        elseif strcmp(file{f},'Ltype')
            datac{n,f} = readpass(folder,info{n,f}(:,3:5));
        else
            datac{n,f} = readpass(folder,info{n,f}(:,end-2:end)); 
        end
       
        %close(findobj(0,'Name', [neuron{n} ' passive']))% 
        %fig = figure('Name',[neuron{n} ' passive'],'NumberTitle','off');
        for e=1:size(info{n,f},1)
            X = linspace(0,1,5000);
            cparam = lsqcurvefit(cfun,cp0,X',datac{n,f}(1:5000,1,e,1),...
                        [-inf,0,-inf],[inf,inf,inf],opts);
            pass(cnt,1:3) = cparam;

            pass(cnt,4) = pass(cnt,1)/range(datac{n,f}(1:5000,2,e,1));
            pass(cnt,5) = pass(cnt,2)/pass(cnt,4);
            pass(cnt,6) = pass(cnt,5)*100;
            group{cnt} = neuron{n};

    %         ax = subplot(2,3,e);
    %         p2 = plot(X,squeeze(datac{n}(:,2,e,1)));hold on
    %         pp = plot(X,cfun(cparam,X),'k');hold on
    %         uistack(p2,'bottom')
    %         ax.YLim = [-140, -70];
    %         ax.Title.String = info{n}{e,7};
            cnt = cnt + 1;
        end    
    end
end


disp('finished loading data')

%% plot properties
disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
figure('Name','scatter','NumberTitle','off');

colors = makecolor;
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];

%boxplots
axs1(1) = axes('Position',[0.05 0.55  0.12 0.4]);
axs1(2) = axes('Position',[0.25 0.55 0.12 0.4]);
axs1(3) = axes('Position',[0.45 0.55 0.12 0.4]);
axs1(4) = axes('Position',[0.65 0.10 0.12 0.4]);
pass2 = pass(:,[4 5 1 2]);
pass3 = prod(pass2,2).^-1;
for n=1:3
    for p=1:4
        if p==3
            Y = pass3(ismember(group,neuron{n}));
        else
            Y = pass2(ismember(group,neuron{n}),p);
        end
        axes(axs1(p))
        xj = histjitter(Y,0.25,15);
        scatter(ones(size(Y))*n+xj,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
axs1(1).YLabel.String = 'Input resistance (M\Omega)';
axs1(2).YLabel.String = 'Capacitance';%'Membrane Area (mm^2)';
axs1(3).YLabel.String = 'Conductance density';
axs1(4).YLabel.String = 'Tau(s)';


disp('Input resistance')
quickanova(pass2(:,1),group,false);
disp('Capacitance')
quickanova(pass2(:,2),group,false);
% disp('Conductance Density (1/gmm2)')
% quickanova(pass3,group,false);
disp('Time constant')
quickanova(pass2(:,4),group,false);

%
thr = cell(1,3);
cnt2 = 1;
logg = '01';
axs1(5) = axes('Position',[0.65 0.55 0.12 0.4]);
%figure
%ax2 = axes;
rheo = cell(1,3);
for n=1:3
    rheo{n} = zeros(1,0);
    thr{n} = nan(150,20);
    cnt = 1;
    %axes(ax2)
    for f=1:3
        if isempty(datac{n,f});continue;end
        for e=1:size(datac{n,f},3)
            I = squeeze(datac{n,f}(1:5000,2,e,:));
            I = I(:)-I(1,1);
            mv = squeeze(datac{n,f}(1:5000,1,e,:));
            mv = mv(:);
            mv = lowpassf(mv,'Fpass',300,'Fstop',700,'Fs',5000,'Apass',1,'Astop',10,'method',{'butter'});
            idx = [diff(mv)>0.25 & [diff(mv,2)>0.25;false] & mv(2:end)>-60;false];
            idx = find(idx);
            
            pulse = I(50:end) - I(1:end-49);
            pup = logg((pulse>0.055 | pulse<-0.055)+1);
            pt = strfind(pup,repelem('01',47));
            sel = arrayfun(@(x) any((x-idx)<10 & (x-idx)>0) & ~any(x-pt<5000 & x-pt>0),idx);
            idx(sel) = [];
            thr{n}(1:length(idx),cnt) = mv(idx);
            if ~isempty(idx);rheo{n} = [rheo{n}, I(idx(1))];end
%             figure;plot(mv);hold on;scatter(idx(:),mv(idx));
%             scatter(ones(size(mv(idx)))*cnt2,mv(idx),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on

            cnt = cnt+1;
            cnt2 = cnt2+1;
        end
    end
    axes(axs1(5))
    %Y = nanmean(thr{n});
    Y = rheo{n}(~isnan(rheo{n}));
    scatter(ones(size(Y))*n,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Y,[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))   
    
end
% axs1(1).YLim
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',neuron);
set(axs1,'XLim',[0.2 3.8]);
axs1(5).YLabel.String = 'Reobase (nA)';


save('pass','datac')

%% Plateau duration

pdur = cell(1,3);
pduri = cell(size(datac));
for n=1:3
    figure("Name",neuron{n})
    cnt = 1;
    pdur{n} = zeros(0,1);
    for k=1:size(datac,2)
        if ~isempty(datac{n,k})
            sump = squeeze(sum(datac{n,k}(:,1,:,:)>-70));
            vals = zeros(size(datac{n,k},3),1);
            for p=1:size(datac{n,k},3)
                idx = find(sump(p,:),1,'last');
                if n==3
                    vals(p) = sump(p,idx) - 25000;
                else
                    vals(p) = sump(p,idx) - 5000;
                end
                ax = subplot(4,9,cnt);
                plot((1:length(datac{n,k}))/5000, squeeze(datac{n,k}(:,1,p,idx)))
                text(10,30,[num2str(round(vals(p)/5000,2)) ' s'])
                ax.XLim = [0 22];
%                 ax.YLim = [-90 50];
%                 ax.YLim = [-2 2];
                if mod(cnt,9)==1
                    ax.YLabel.String = 'mV';
                end

                if cnt>27
                    ax.XLabel.String = 'Time (s)';
                end
                cnt = cnt + 1;
            end
            pdur{n} = [pdur{n}; squeeze(vals)/5000];
        end
    end
    pdur{n}(pdur{n}<0) = 0;
end

%%
figure;
subplot(1,2,1)
durt = 0.5;
pie([sum(pdur{1}>durt), sum(pdur{1}<=durt)]/length(pdur{1}));
legend(["plateau","non-plateau"],'Location','south')
subplot(1,2,2)
pie([sum(pdur{2}>durt), sum(pdur{2}<=durt)]/length(pdur{2}));
legend(["plateau","non-plateau"],'Location','south')



