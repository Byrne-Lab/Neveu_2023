% 
% folder = 'D:\Data\2019\';
% 
% neuron = ["B51","B64","B8"];
% info=cell(1,3);
% data = cell(1,3);
% datac = cell(1,3);
% command = cell(1,3);
% logg = '01';
% for n=1:length(neuron)
%     [~,~,info{n}] = xlsread([folder 'KCa_new_4mTEA_NiCd.xlsx'],neuron{n});
%     info{n} = info{n}(3:end,2:end);
%     info{n} = cellfun(@num2str,info{n},'UniformOutput',false);
%     if ~isempty(info{n})
%         [data{n},si] = readsdata(folder,info{n}(:,1:2),1:1.1e4,4.5e3);
%         command{n} = round(squeeze(data{n}(5000,2,:,3,:))/10)'*10;
%         data{n}(:,[1 3],:,3,:) = lowpassf(data{n}(:,[1 3],:,3,:),'Fpass',700,'Fstop',900,'Fs',1e6/si);% 500 700
%         
%         for e=1:size(info{n},1)
%             [datagf,sig] = abfload([folder info{n}{e,3} '.abf'],'verbose',false);
%             datagf = datagf(round(str2double(info{n}{e,4})*1e3/sig):round(str2double(info{n}{e,5})*1e3/sig),:);
%             pulse = datagf(50:end,4) - datagf(1:end-49,4);
%             pup = logg((pulse>0.055)+1);
%             pt = strfind(pup,repelem('01',47));
%             pdn = logg((pulse<-0.055)+1);
%             pe = strfind(pdn,repelem('01',47));
%             long = arrayfun(@(x) any(abs(pt-x-25000)<100),pe);       
%             pe = pe(long);
%             datac{n}(:,:,e) = datagf(pe(1)+50:pe(1)+5000,:);
%         end
%     end
% end
% 
% 
% disp('>>>> finished reading data <<<<<')
% 
% save('KCa_NiCd','data','command','datac','info','si')

%% current traces

neuron = ["B51","B64","B8"];

load('KCa_NiCd','data','command','datac','info','si');

for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtractions']))% 
    fig = figure('Name',[neuron{n} ' subtractions'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(2,4,e);
%         plot(400:0.1:1000,squeeze(mean(data{n}(4000:10000,1,:,3,:),5))) 
        Y = squeeze(data{n}(4000:10000,1,:,3,e));
        plot(400:0.1:1000,Y)  
%         ax.YLim = [-40, 250];
%         ax.YLim = [-50, 50];
        ax.YLim = [-50, 180];
        ax.XLabel.String = 'Time (ms';
        ax.Title.String = info{n}{e,1};%['Pre: ', num2str(info{e,3}), ', Post: ', num2str(info{e,4}),', Pre: ', num2str(info{e,5})];
    end
end


%% current traces

colors = makecolor;
dcolors = makecolor(-0.3);
bcolors = makecolor(0.5);

close(findobj(0,'Name', 'subtractions'))% 
fig = figure('Name','subtractions','NumberTitle','off');
for n=1:length(neuron)
            ax = subplot(1,3,n);
    for t=1:9

        Y = squeeze(data{n}(4000:10000,1,t,3,:));
        X = 1:size(Y,1);
        mY = mean(Y,2);
        
        fx = [X fliplr(X)];
        err = nanstd(Y,[],2)/sqrt(size(Y,2));
        err = [-err';err'] + nanmean(Y');
        fy = [err(1,:) fliplr(err(2,:))];
        fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(n,:));hold on
        
        plot(mY,'Color',dcolors(n,:));hold on 
        
    %         ax.YLim = [-40, 180];
    %     ax.YLim = [-90, 150];
    %     ax.Title.String = info{n}{e,1};%['Pre: ', num2str(info{e,3}), ', Post: ', num2str(info{e,4}),', Pre: ', num2str(info{e,5})];
    end
end

%% passive properties
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
pass = cell(1,3);
cfun = @(p,x) p(1).*exp(-x./p(2)) + p(3);
cp0 = [20,0.2,-80];
for n=1:length(neuron)
    %close(findobj(0,'Name', [neuron{n} ' passive']))% 
    %fig = figure('Name',[neuron{n} ' passive'],'NumberTitle','off');
    for e=1:size(info{n},1)
        X = linspace(0,1,size(datac{n},1));
        cparam = lsqcurvefit(cfun,cp0,X',datac{n}(:,2,e,1),...
                    [-inf,0,-inf],[inf,inf,inf],opts);
        pass{n}(e,1:3) = cparam;

        pass{n}(e,4) = pass{n}(e,1)/range(datac{n}(:,3,e,1));
        pass{n}(e,5) = pass{n}(e,2)/pass{n}(e,4);
        pass{n}(e,6) = pass{n}(e,5)*100;

%         ax = subplot(2,3,e);
%         p2 = plot(X,squeeze(datac{n}(:,2,e,1)));hold on
%         pp = plot(X,cfun(cparam,X),'k');hold on
%         uistack(p2,'bottom')
%         ax.YLim = [-140, -70];
%         ax.Title.String = info{n}{e,7};
    end
end

%% fits
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
io = 4440;
fW = 4600:9380;
delta = (fW(1) - io)*si/1e6;
X = linspace(0,range(fW)*si/1e6,length(fW)) + delta;
KAp = cell(1,3);
KAa = cell(1,3);

disp('fits')
for n=1:length(neuron)
    for e=1:size(info{n},1)
        disp([neuron{n},'  ',info{n}{e,1}])
        for t=1:size(data{n},3)
            if ~isnan(data{n}(1,1,t,3,e))
                yy = data{n}(fW,1,t,3,e);

                fun = @(p,x) p(1).*(1-p(3).*(1-exp(-x./p(2))));% for single exponential
                p0 = [1     0.2     0];
                lb = [0     0.05    0];
                ub = [inf   2       1];
                KAp{n}(e,t,:) = lsqcurvefit(fun,p0,X(1:length(fW)),yy',lb,ub,opts);

            else
                KAp{n}(e,t,:) = nan;
            end
        end
    end
    KAa{n}(KAa{n}==Inf) = 0;
end



KAa = cell(1,3);
KAg = cell(1,3);
for n=1:3
    gg = nan(size(info{n},1),size(data{n},3));
    for e=1:size(info{n},1)
        for t=1:size(data{n},3)
            gg(e,t) = KAp{n}(e,t,1) ./(command{n}(e,t) + 70);
        end
    end
    gg(gg>10) = 0;
    m = max(gg,[] ,2);
    KAg{n} = m;
    m = repmat(m,1,length(command{n}));
    KAa{n} = gg./m;%KAp{n}(:,:,1)./(m.*(command{n} + 70));
end



%% check fits

expp = 1:5;

bln = 400;
fW3 = fW(1)-bln:fW(end);%6304;
X = linspace(0,range(fW)*si/1e6,length(fW));
X2 = [fliplr(-X(1:bln)) , X];


colors = makecolor;
dcolors = makecolor(-0.3);
for n=1:3
    close(findobj(0,'Name', [neuron{n} 'fits']))% 
    figure('Name',[neuron{n} 'fits'],'NumberTitle','off','Position',[50+(n-1)*560 132 560 770]);
    cnt = 1;
    for e=1:size(info{n},1)
        if e>size(KAp{n},1);continue;end
        axp = subplot(2,ceil(length(expp)/2),cnt);
        cnt=cnt+1;
        plt = gobjects(0,1);
        for t=4:size(data{n},3)
            plt(end+1) = plot(X2,data{n}(fW3,1,t,3,e),'Color',colors(t-3,:));hold on
            plot(X,fun(KAp{n}(e,t,:),X),'Color',dcolors(t-3,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-10 110];
    end
    legend(plt,num2str(command{n}(1,4:end)'))
end


%% scatters
% ---------------------------------------------------------------------------------
disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
figure('Name','scatter','NumberTitle','off','Position',[1000 100 1000 800]);
names = ["B51","B64","B8"];

% activation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
axs0 = axes('Position',[0.05 0.55 0.3 0.4]);
colors = makecolor;
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];

bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2)));
bp0 = [-10, 1];
KAap = cell(1,3);
aplt = gobjects(1,3);
for n=1:3
    iqr = prctile(KAa{n},[25 75]);
    err = nanstd(KAa{n})/sqrt(size(KAa{n},1));
    err = [-err;err] + nanmean(KAa{n});

    Xx = command{n}(1,:)+n-2;
    scatter(Xx,nanmean(KAa{n}),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    for e=1:size(KAp{n},1)
        isn = ~isnan(KAa{n}(e,:));
        X = command{n}(e,isn);
        Y = KAa{n}(e,isn);
        KAap{n}(e,:) = lsqcurvefit(bfun,bp0,X(2:end),Y(2:end),[-inf,-inf],[inf,inf],opts);
    end
    aplt(n) = plot(-80:20,bfun(mean(KAap{n}),-80:20),'Color',colorss(n,:));hold on
end
legend(aplt,names,'Location','northwest')
axs0.Title.String = 'Activation Curve';
axs0.YLabel.String = 'Activation';
axs0.XLabel.String = 'Command potential (mV)';
axs0.XLim = [-80,20];
axs0.YLim = [0 1];
axs0.TickDir = 'out';

%boxplots
axs1(1) = axes('Position',[0.4 0.55  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.55 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.55  0.12 0.4]);
Yah = zeros(0,2);
Ygh = cell(0,1);
YaA = zeros(0,2);
for n=1:3
    Yah = [Yah;KAap{n}]; %#ok<AGROW>
    Ygh = [Ygh;repmat({names{n}},size(KAap{n},1),1)];
    for p=1:2
        axes(axs1(p))
        scatter(ones(size(KAap{n},1),1)*n,KAap{n}(:,p),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(KAap{n}(:,p),[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
    axes(axs1(3))
    Ys = KAp{n}(:,end,1)/(command{n}(1,end) +70);%./pass{n}(:,6);
    YaA = [YaA; Ys];
    scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Ys,[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));
end 
axs1(1).YLabel.String = 'Half activation (mV)';
axs1(2).YLabel.String = 'Slope';
axs1(3).YLabel.String = 'Maximum conductance (g)';%'Current Density (nA/mm^2)';
axs1(2).YLim(1) = 0;
axs1(3).YLim(1) = 0;
set(axs1,'TickDir','out');
% axs1(3).YLim(2) = 3.5;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

disp('Half activation')
[cta{1},paa(1)] = quickanova(Yah(:,1),Ygh);
makesig(cta{1},paa(1),axs1(1));
disp('Slope')
[cta{1},paa(1)] = quickanova(Yah(:,2),Ygh);
makesig(cta{1},paa(1),axs1(2));
disp('Current Density')
[cta{1},paa(1)] = quickanova(YaA,Ygh);
makesig(cta{1},paa(1),axs1(3))








% inactivation
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
axs2(1) = axes('Position',[0.05 0.05 0.3 0.4]);

bfun = @(p,x) (1 - p(3))./(1+exp((p(1) - x)/-5)) + p(3);
bp0 = [-10, -5,0];
KAiap = cell(1,3);%KAp{n}(e,t,:,s)
for n=1:length(neuron)
    idx = 6:9;
    Ya = KAp{n}(:,idx,3);
    Xx = command{n}(1,idx);
    iqr = prctile(Ya,[25 75]);
    err = nanstd(Ya)/sqrt(size(KAp{n},1));
    err = [-err;err] + 1-nanmean(Ya);

    
    plot(repmat(Xx+n-2,2,1),err,'Color',colorss(n,:));hold on
    scatter(Xx+n-2,1-nanmean(Ya),20,'d','MarkerFaceColor','w','MarkerEdgeColor',colorss(n,:));hold on
    for e=1:size(Ya,1)
        isn = ~isnan(Ya(e,:));
        Xc = [-60, Xx(isn)];
        Yi = [1, 1-Ya(e,isn)];
        KAiap{n}(e,:) = lsqcurvefit(bfun,bp0,Xc,Yi,[-inf,-inf,0],[inf,inf,1],opts);
    end
    iaplt(n) = plot(-80:20,bfun(mean(KAiap{n}),-80:20),'Color',colorss(n,:),'LineStyle','--');hold on
end
legend(iaplt,names,'Location','northwest')
axs2(1).Title.String = 'Inactivation Curve';
axs2(1).YLabel.String = 'B';
axs2(1).XLabel.String = 'Command potential (mV)';
axs2(1).YLim = [0 1];
axs2(1).XLim = [-80,21];


%boxplots
axs1(1) = axes('Position',[0.4 0.05  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.05 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.05  0.12 0.4]);
Yiah = zeros(0,2);
Yigh = cell(0,1);
YiaA = zeros(0,2);%KAp{n}(e,t,:,s)
for n=1:length(neuron)
    Yiah = [Yiah;KAiap{n}]; %#ok<AGROW>
    Yigh = [Yigh;repmat({names{n}},size(KAiap{n},1),1)];
    for p=1:3
        axes(axs1(p))
        scatter(ones(size(KAiap{n},1),1)*n,KAiap{n}(:,p),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(KAiap{n}(:,p),[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
axs1(1).YLabel.String = 'Half inactivation (mV)';
axs1(2).YLabel.String = 'Slope';
axs1(3).YLabel.String = 'Maximum Inactivation';
axs1(3).YLim = [0 1];
axs1(2).YLim(2) = 0;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

disp('==================Inactivation================')
disp('Half inactivation')
[cta{1},paa(1)] = quickanova(Yiah(:,1),Yigh);
makesig(cta{1},paa(1),axs1(1));
disp('Inactivation Slope')
[cta{1},paa(1)] = quickanova(Yiah(:,2),Yigh);
makesig(cta{1},paa(1),axs1(2));
disp('Maximum Inactivation')
[cta{1},paa(1)] = quickanova(Yiah(:,3),Yigh);
makesig(cta{1},paa(1),axs1(3));




%% scatter time constant
close(findobj(0,'Name', 'scatter time constant'))% 
figure('Name','scatter time constant','NumberTitle','off');

% fit inactivation time constant


atp = cell(1,3);
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
% mint = 0.001;
for n=1:length(neuron)
    for e=1:size(KAp{n},1)
        Yd = KAp{n}(e,idx,2);
        xt = command{n}(e,idx);
        tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(4))) + p(2);
        tp0 = [0.1, 0, -20, 1];%max(Yd)
    end
end

axs22 = axes('Position',[0.05 0.55 0.4 0.4]);
atplt = gobjects(1,3);
for n=1:length(neuron)
    Ym = mean(KAp{n}(:,idx,2));
    err = nanstd(KAp{n}(:,idx,2))/sqrt(size(KAp{n}(:,idx,2),1));
    Ye = [Ym + err;Ym - err];
    plot(command{n}(1:2,idx)+n-2,Ye,'Color',colorss(n,:));hold on
    scatter(command{n}(1,idx)+n-2,Ym,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    switch n
        case 1
            atplt(n) = plot(-80:30,tfun([0.8 0.08 18 2 20],-80:30),'Color',colorss(n,:));hold on
        case 2
            atplt(n) = plot(-80:30,tfun([1.1 0.2 -13 2 -6],-80:30),'Color',colorss(n,:));hold on
        case 3
            atplt(n) = plot(-80:30,tfun([1.4 0.2 -2  1.5  3.5],-80:30),'Color',colorss(n,:));hold on
    end
end
% legend(names)
axs22.Title.String = 'Inactivation time constant';
% axs22.YLim(1) = 0;
axs22.YLabel.String = 'Time constant (s)';
axs22.XLim(1) = -60;



% fit inactivation time constant scatter


axs23 = axes('Position',[0.55 0.55 0.4 0.4]);
for n=3%1:length(neuron)
    for e=1:size(KAp{n},1)
        Y = KAp{n}(e,idx,2);
        scatter(command{n}(1,idx),Y,20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    plot(command{n}(1,idx),mean(KAp{n}(:,idx,2)))
end




