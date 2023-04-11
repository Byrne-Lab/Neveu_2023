folder = '\\uthsch-nas.uthouston.edu\ms_nba\ByrneLab\cneveu\Data\2019\';
if strcmp(getenv('COMPUTERNAME'),'MSNBA160674') || strcmp(getenv('COMPUTERNAME'),'MSNBA-C000045')
%     folder = 'C:\Users\cneveu\Desktop\Data\2019\';
%     folder = 'D:\Data\2019\';
%     folder = 'E:\external_drive\Data\2019\';
%     folder = 'E:\Data\2019\';
    folder = 'C:\Users\cneveu\My Drive\Data\2019\';
end
% folder = 'E:\external_drive\Data\2019\';

neuron = ["B51","B64","B8"];
info=cell(1,3);
data = cell(1,3);
datai = cell(1,3);
datac = cell(1,3);
datar = cell(1,3);
vrange = [-30 , 20];
command = cell(1,3);
commandi = cell(1,3);
xidx = [1160 4300];
for n=1:length(neuron)
    [~,~,info{n}] = xlsread([folder 'Delayed_K.xlsx'],neuron{n});
    info{n} = info{n}(3:end,2:end);
    info{n} = cellfun(@num2str,info{n},'UniformOutput',false);

    [data{n},si] = readsdata(folder,info{n}(:,[1 5]),1:2e4,4.5e3,10);
    data{n}(:,[1 3],:,1,:) = lowpassf(data{n}(:,[1 3],:,1,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
    
    
    command{n} = round(squeeze(data{n}(5000,2,:,1,:))'/10)*10;
    if any(command{n}(:) ~= repelem(command{n}(1,:)',size(command{n},1))); error('pulses different'); end
    sel = command{n}(1,:)>=vrange(1) & command{n}(1,:)<=vrange(2);
    command{n} = command{n}(:,sel);
    data{n}(:,:,~sel,:,:) = [];
   
    sz = size(data{n});
%     leakd = nan([1,1,1,2,sz(5),2]);
%     for i=1:length(xidx)
%         leakd(:,:,:,:,:,i) = mean(data{n}(xidx(i)-20:xidx(i)+10,1,:,1:2,:),[1 3]);
%     end
%     leakm = diff(leakd,1,6)/40;
%     sub =  - repmat(leakm,[sz(1),1,sz(3),1,1]).*(data{n}(:,2,:,1:2,:)+80)  - repmat(leakd(:,:,:,:,:,1),[sz(1),1,sz(3),1,1]);
%     data{n}(:,1,:,1:2,:) = data{n}(:,1,:,1:2,:) - repmat(leakm,[sz(1),1,sz(3),1,1]).*(data{n}(:,2,:,1:2,:)+80)...
%                                                 - repmat(leakd(:,:,:,:,:,1),[sz(1),1,sz(3),1,1]);
    
%     % recovery
%     for e=1:size(info{n},1)
%         [datar{n}(:,:,:,1,e),sir] = abfload([folder info{n}{e,3} '.abf'],'verbose',false);
%         [datar{n}(:,:,:,2,e),~ ] = abfload([folder info{n}{e,6} '.abf'],'verbose',false);
%     end
%     datar{n}(:,:,:,3,:) = datar{n}(:,:,:,1,:) - datar{n}(:,:,:,2,:);
%     datar{n}(:,2,:,3,:) = datar{n}(:,2,:,1,:);
%     datar{n}(:,[1 3],:,3,:) = lowpassf(datar{n}(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
%     datar{n}(10001:end,:,:,:,:) = [];
    
    dta = readpass(folder,info{n}(:,9:11)); 
    datac{n} = dta(:,:,:,1);
end

data{3}(5761,[1 3],1,:,1) = data{3}(5760,[1 3],1,:,1);
disp('>>>> finished reading data <<<<<')

%% passive properties
disp('passive')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
pass = cell(1,3);
cfun = @(p,x) p(1).*exp(-x./p(2)) + p(3);
cp0 = [20,0.2,-80];
X = linspace(0,1,5000);
for n=1:length(neuron)
    %close(findobj(0,'Name', [neuron{n} ' passive']))% 
    %fig = figure('Name',[neuron{n} ' passive'],'NumberTitle','off');
    for e=1:size(info{n},1)
        cparam = lsqcurvefit(cfun,cp0,X',datac{n}(1:5000,1,e,1),...
                    [-inf,0,-inf],[inf,inf,inf],opts);
        pass{n}(e,1:3) = cparam;

        pass{n}(e,4) = pass{n}(e,1)/range(datac{n}(:,2,e,1));
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



%% Activation current traces
s=1;
colors = makecolor;
dcolors = makecolor(-0.3);
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtractions']))% 
    fig = figure('Name',[neuron{n} ' subtractions'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(3,3,e);
        for t=1:size(data{n},3)
            plot(squeeze(data{n}(:,1,t,s,e)),'Color',colors(t,:));hold on
        end
        ax.YLim = [-30, 300];
        ax.Title.String = info{n}{e,1};
    end
end

%% inact fits
disp('run fits')
opts1 = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
fWi = 14500:15250;
X = linspace(0,length(fWi),length(fWi));
DKpi = cell(1,3);
funi = @(p,x) p(1).*exp(-x./p(2)) + p(3);
measin = cell(1,3);
back = -190;
ffig = gobjects(1,3);
sfig = gobjects(1,3);

opts2 = optimset('Display','off','Algorithm','levenberg-marquardt');
iifun = @(p,x) (p(1) - p(2))./(1+exp((p(3) - x)/p(4))) + p(2);
ip0 = [200 50 -20 -5];
ifits = cell(1,3);
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} 'inact fits']))% 
    ffig(n) = figure('Name',[neuron{n} 'inact fits'],'NumberTitle','off');
    close(findobj(0,'Name', [neuron{n} 'scatter inact fits']))% 
    sfig(n) = figure('Name',[neuron{n} 'scatter inact fits'],'NumberTitle','off');
    for e=1:size(info{n},1)
        figure(ffig(n))
        axp = subplot(2,3,e);
        plt = gobjects(1,size(data{n},3));
        disp([neuron{n},'  ',info{n}{e,1}])
        for t=1:size(data{n},3)
            yy = data{n}(fWi,1,t,s,e);
            p0 = [ max(yy) , 0.1  , data{n}(fWi(end),1,t,s,e)];
%             DKpi{n}(e,t,:,s) = lsqcurvefit(funi,p0,X',yy,[0,300,0],[Inf,Inf,Inf],opts1);
%             measin{n}(t,s,e) = funi(DKpi{n}(e,t,:,s),back);
            measin{n}(t,s,e) = yy(1);
            fWi2 = 1:16000;
            plt(t) = plot(fWi2,data{n}(fWi2,1,t,s,e),'Color',colors(t,:));hold on
%             plot(X+back+fWi(1)-fWi2(1),funi(DKpi{n}(e,t,:,s),X+back),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        scatter(fWi(1),data{n}(fWi(1),1,t,s,e)+10,'+')
        
        ifits{n}(e,:) = lsqcurvefit(iifun,ip0,command{n}(e,:),measin{n}(:,s,e)',[-inf,-inf,-inf,-inf],[Inf,Inf,Inf,inf],opts2);
        figure(sfig(n))
        axxs = subplot(2,3,e);
        scatter(command{n}(e,:),measin{n}(:,s,e)/ifits{n}(e,1),'MarkerEdgeColor',colorss(n,:));hold on
        plot(-40:20,iifun(ifits{n}(e,:),-40:20)/ifits{n}(e,1),'k');hold on
        axxs.YLim(1) = 0;
%         measin{n}(:,s,e) = 1 - measin{n}(:,s,e)/measin{n}(1,s,e);
        measin{n}(:,s,e) = 1 - iifun(ifits{n}(e,:),command{n}(e,:))/ifits{n}(e,1);
%         axp.Title.String = info{n}{e,1};
%         axp.YLim = [0 600];
        axp.XLim = [14300 15200];
    end
    legend(plt,num2str(command{n}(1,:)'))
end

%% fits
disp('run fits')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
% measin = cell(1,3);
fW = 4320:14300;
X = linspace(0,length(fW)/1e4,length(fW));
DKp = cell(1,3);
ofst = 60;
for n=1:length(neuron)
%     measin{n} = squeeze(mean(data{n}(14500:14520,1,:,:,:)));% 14200:14300
%     init = mean(data{n}(14500:14520,1,1,:,:));
%     init = repmat(init,[1,1,size(data{n},3),1,1]);
%     measin{n} = ones(size(measin{n})) - measin{n}./squeeze(init);%measin{n}(t,s,e)
%     measin{n}(measin{n}<0) = 0;
    for e=1:size(info{n},1)
        disp([neuron{n},'  ',info{n}{e,1}])
%         for s=1:3
            for t=1:size(data{n},3)
                yy = data{n}(fW(1+ofst:end),1,t,s,e);
%                 fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - measin{n}(t,s,e).*(1-exp(-x./p(4))));
%                 fun = @(p,x) p(1).*(1 - exp(-x./0.01)).^2.*(1 - measin{n}(t,s,e).*(1-exp(-x./(0.01+p(4)))));
                fun = @(p,x) p(1).*(1 - exp(-x./0.01)).*(1 - measin{n}(t,s,e).*(1-exp(-x./(0.01+p(4)))));
                p0 = [ prctile(yy,99)   , 0.01, measin{n}(t,s,e),  2];
                lb = [prctile(yy,99)    , 0   ,     0,            0.15 ];
                ub = [prctile(yy,99)*1.5, 0.1 ,     1,            Inf];
                DKp{n}(e,t,:,s) = lsqcurvefit(fun,p0 ,X(1+ofst:end)',yy,lb,ub,opts);
                DKp{n}(e,t,4,s) = sum(DKp{n}(e,t,[2 4],s));
                DKp{n}(e,t,3,s) = measin{n}(t,s,e);
            end
%         end
    end
end
%

fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
KAa = cell(1,3);
maxg = cell(1,3);
for n=1:3
    gg = squeeze(DKp{n}(:,:,1,:)) ./repmat((command{n} +85),[1,1]);
    gg(gg>10) = 0;
    m = max(gg,[] ,2);
    maxg{n} = m;
    m = repmat(m,1,length(command{n}));
    KAa{n} = gg./m;%KAp{n}(:,:,1)./(m.*(command{n} + 70));
end

%% check fits

bln = 400;
fW = 4320:14300;
fW3 = fW(1)-bln:fW(end);%6304;
X = linspace(0,range(fW)*si/1e6,length(fW));
X2 = [fliplr(-X(1:bln)) , X];
exxp = [1 2 5];

red = 4;
fW3 = fW3(1:red:end);
X2 = X2(1:red:end);
X = X(1:red:end);

% X = linspace(0,length(fW)/1e4,length(fW));
% fun = @(p,x) p(1).*(1 - exp(-x./p(2))).^2.*(1 - p(3).*(1 - exp(-x./p(4))));
fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - p(3).*(1 - exp(-x./p(4))));
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} 'fits']))% 
    figure('Name',[neuron{n} 'fits'],'NumberTitle','off');%,'Position',[(n-1)*600 100 600 900]);
    cnt = 1;
    for e=1:size(info{n},1)
        axp = subplot(2,3,cnt);
%         axp = subplot(1,1,cnt);
        plt = gobjects(1,size(data{n},3));
        for t=1:size(data{n},3)
            plt(t) = plot(X2,data{n}(fW3,1,t,s,e),'Color',colors(t,:));hold on
            plot(X,fun(DKp{n}(e,t,:,s),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-100 350];
        cnt = cnt+1;
    end
    legend(plt,num2str(command{n}(1,:)'))
end


%% scatters
disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
figure('Name','scatter','NumberTitle','off');%,'Position',[1000 500 1000 800]);
names = ["B51","B64","B8"];

% activation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
axs0 = axes('Position',[0.05 0.55 0.3 0.4]);
colors = makecolor;

bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2))).^2;
bp0 = [-10, 1];
KAap = cell(1,3);
for n=1:length(neuron)
%     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    Yy = KAa{n}(:,:,s);
%     Yy = Yy - repmat(min(KAa{n}(:,:,s),[],2),1,size(Yy,2));
    iqr = prctile(Yy,[25 75]);
    err = nanstd(Yy)/sqrt(size(KAa{n},1));
    err = [-err;err] + nanmean(Yy);

    Xx = command{n}(1,:)+n-2;
%     scatter(Xx,nanmean(Yy),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
%     plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    for e=1:size(DKp{n},1)
        isn = ~isnan(KAa{n}(e,:,s));
        Xc = command{n}(e,isn);
        Y = Yy(e,isn);
        KAap{n}(e,:) = lsqcurvefit(bfun,bp0,Xc,Y,[-inf,-inf],[inf,inf],opts);
%         plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%         scatter(command{n}(e,:)+n,KAa{n}(e,:),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    scatter(Xx,nanmean(Yy),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    aplt(n) = plot(-80:20,bfun(mean(KAap{n}),-80:20),'Color',colorss(n,:));hold on
end
axs0.Title.String = 'Activation Curve';
axs0.YLabel.String = 'Activation';
axs0.XLabel.String = 'Command potential (mV)';
axs0.XLim = [-80,21];


%boxplots
axs1(1) = axes('Position',[0.4 0.55  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.55 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.55  0.12 0.4]);
Yah = zeros(0,2);
Ygh = cell(0,1);
YaA = zeros(0,2);%DKp{n}(e,t,:,s)
for n=1:length(neuron)
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
%     YaA = [YaA; DKp{n}(:,end,1)./pass{n}(:,6)];
%     scatter(ones(size(DKp{n},1),1)*n, DKp{n}(:,end,1,s)./pass{n}(:,6),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
%     Yp = prctile(DKp{n}(:,end,1)./pass{n}(:,6),[25 50 75]);
%     rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
%     plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));
    Ys = DKp{n}(:,end,1)/(command{n}(1,end) +85);%./pass{n}(:,6);
    YaA = [YaA; Ys];
    scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Ys,[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));
end 
axs1(1).YLabel.String = 'Half activation (mV)';
axs1(2).YLabel.String = 'Slope';
axs1(3).YLabel.String = 'Current Density (nA/mm^2)';
axs1(3).YLim(1) = 0;
axs1(2).YLim(1) = 0;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

disp('Half activation')
[cta{1},paa(1)] = quickanova(Yah(:,1),Ygh);
makesig(cta{1},paa(1),axs1(3));
disp('Slope')
[cta{1},paa(1)] = quickanova(Yah(:,2),Ygh);
makesig(cta{1},paa(1),axs1(3));
disp('Current Density')
[cta{1},paa(1)] = quickanova(YaA,Ygh);
makesig(cta{1},paa(1),axs1(3));



% inactivation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
% axs2(1) = axes('Position',[0.05 0.05 0.3 0.4]);
axes(axs0)

bfun = @(p,x) (1 - p(3))./(1+exp((p(1) - x)/p(2))) + p(3);
bp0 = [-10, -1,0];
KAiap = cell(1,3);%DKp{n}(e,t,:,s)
for n=1:length(neuron)
%     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    iqr = prctile(DKp{n}(:,:,3,s),[25 75]);
    err = nanstd(DKp{n}(:,:,3,s))/sqrt(size(DKp{n},1));
    err = [-err;err] + 1-nanmean(DKp{n}(:,:,3,s));

    Xx = command{n}(1,:)+n-2;
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    scatter(Xx,1-nanmean(DKp{n}(:,:,3,s)),20,'d','MarkerFaceColor','w','MarkerEdgeColor',colorss(n,:));hold on
    for e=1:size(DKp{n},1)
        isn = ~isnan(DKp{n}(e,:,3,s));
        Xc = command{n}(e,isn);
        Y = 1-DKp{n}(e,isn,3,s);
        KAiap{n}(e,:) = lsqcurvefit(bfun,bp0,Xc,Y,[-inf,-inf,-inf],[inf,inf,inf],opts);
%         plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%         scatter(command{n}(e,:)+n,KAa{n}(e,:),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    iaplt(n) = plot(-80:20,bfun(mean(KAiap{n}),-80:20),'Color',colorss(n,:),'LineStyle','--');hold on
end
legend(aplt,names,'Location','southeast')
% legend(iaplt,names,'Location','southeast')
% axs2(1).Title.String = 'Inactivation Curve';
% axs2(1).YLabel.String = 'Activation';
% axs2(1).XLabel.String = 'Command potential (mV)';
% axs2(1).YLim = [0 1];
% axs2(1).XLim = [-60,21];


%boxplots
axs1(1) = axes('Position',[0.4 0.05  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.05 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.05  0.12 0.4]);
Yiah = zeros(0,2);
Yigh = cell(0,1);
YiaA = zeros(0,2);%DKp{n}(e,t,:,s)
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


%% Relaxation from Innactivation

recovfit = false;

opts = optimset('Display','off','Algorithm','levenberg-marquardt');
measp = cell(1,3);
times = cell(1,3);
for n=1:length(neuron)
    if recovfit
        close(findobj(0,'Name', ['recoveryfit' neuron{n}]))
        figure('Name',['recoveryfit' neuron{n}],'NumberTitle','off'); 
    end
    for e=1:size(info{n},1)
        [datas,si,head] = abfload([folder info{n}{e,4} '.abf']);
        if e==1
            times{n} = zeros(size(info{n},1),3*size(datas,3));
            measp{n} = times{n};
        end
        strt = round(head.sweepStartInPts - head.sweepStartInPts(1) + 1);
        cnt = 1;
        for v=1:length(strt)
            logstr = char(double(datas(:,2,v)'>0)+'0');      
            tms = strfind(logstr,'01');
            times{n}(e,(1:3)+(v-1)*3) = tms+strt(v);
            for t=1:length(tms)
                y = datas((0:1e4)+tms(t),3,v);
                strtt = 1e3;
                yf = y(strtt:end);
                pi0 = [range(yf) 5000];
                C = mean(y(end-100:end));
                efun = @(p,x) p(1).*exp(-x./p(2))+C;
                pfip = lsqcurvefit(efun,pi0,(1:length(yf))',yf,-Inf(1,2),Inf(1,2),opts);
                measp{n}(e,t+(v-1)*3) = efun(pfip,-strtt); 
                if recovfit
                    ax(cnt) = subplot(3,4,cnt);
                    plot(y);hold on
                    plot(efun(pfip,-strtt:length(yf)),'k');hold on
                    cnt = cnt+1;
                end
            end
        end
    end
    if recovfit; set(ax,'YLim',[0 200]);end
    times{n} = [zeros(size(times{n},1),1)+60 ,  diff(times{n},1,2)*si/1e6 - 1]; 
    measp{n} = measp{n}./repmat(max(measp{n},[],2),1,12);
end
%%

close(findobj(0,'Name', 'scatter time constant'))% 
figure('Name','scatter time constant','NumberTitle','off');

axs211 = axes('Position',[0.05 0.55 0.4 0.4]);
prec = zeros(3,1);
leg = strings(3,1);
rplt = gobjects(1,3);
for n=1:length(neuron)
    
    subplot(axs211)
    
    %scatter(times(:),measp(:),'kd','MarkerFaceColor','k');hold on
    scatter(times{n}(:),measp{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
%     fexp = @(p,x) 1-0.45.*exp(-x/p(1)) - 0.15.*exp(-x/p(2));
%     prec(n,:) = lsqcurvefit(fexp,[1.1 30],times(:),measp(:),-Inf(1,2),Inf(1,2),opts);
    fexp = @(p,x) 1-(1-mean(KAiap{n}(:,3))).*exp(x/p(1));
    prec(n) = lsqcurvefit(fexp,-4,times{n}(:),measp{n}(:),-Inf,Inf,opts);
    rplt(n) = plot(0:0.1:60,fexp(prec(n),0:0.1:60),'Color',colors(n,:));hold on
    leg(n) = join([neuron{n}," ",num2str(prec(n),3)]);
end
axs211 = gca;
axs211.YLim = [0.4 1];
axs211.XLim(1) = 0;
axs211.XLabel.String = 'Time (s)';
axs211.YLabel.String = 'Activation';
legend(rplt,leg)


atp = cell(1,3);
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
tfun = @(p,x) (0.1 - p(3))./(1+exp((p(1) - x)/p(2))) + p(3);
tp0 = [-15 -4 0.005];
for n=1:length(neuron)
    for e=1:size(DKp{n},1)
        Ya = DKp{n}(e,:,2);
        atp{n}(e,:) = lsqcurvefit(tfun,tp0,command{n}(e,:),Ya,[-inf,-inf,-inf],[inf,inf,inf],opts);
    end
end

% time constant scatter
axs21 = axes('Position',[0.05 0.05 0.4 0.4]);
atplt = gobjects(1,3);
for n=1:length(neuron)
    Ym = mean(DKp{n}(:,:,2));
    err = nanstd(DKp{n}(:,:,2))/sqrt(size(DKp{n}(:,:,2),1));
    Ye = [Ym + err;Ym - err];
    plot(command{n}(1:2,:)+n-2,Ye,'Color',colorss(n,:));hold on
    scatter(command{n}(1,:)+n-2,Ym,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    atplt(n) = plot(-40:0.1:20,tfun(mean(atp{n}),-40:0.1:20),'Color',colorss(n,:));hold on
end
legend(atplt,names)
axs21.Title.String = 'Activation time constant';
axs21.YLim = [0 0.11];
axs21.YLabel.String = 'Time constant (s)';


%boxplots
axs43(1) = axes('Position',[0.5 0.05  0.09 0.4]);
axs43(2) = axes('Position',[0.62 0.05  0.09 0.4]);
axs43(3) = axes('Position',[0.74 0.05  0.09 0.4]);
Yiah = zeros(0,2);
Yigh = cell(0,1);
YiaA = zeros(0,3);%DKp{n}(e,t,:,s)
for n=1:length(neuron)
    Yiah = [Yiah;atp{n}]; %#ok<AGROW>
    Yigh = [Yigh;repmat({neuron{n}},size(atp{n},1),1)]; %#ok<AGROW>
    for p=1:size(atp{n},2)
        axes(axs43(p))
        Y = atp{n}(:,p);
        scatter(ones(size(atp{n},1),1)*n,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        text(ones(size(atp{n},1),1)*n+0.1,Y,string((1:length(Y))'))
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
set(axs43,'XTick',1:3);
set(axs43,'XTickLabel',neuron);
set(axs43,'XLim',[0.2 3.8]);
% set(axs43([1 2]),'YLim',[0 inf]);

ylab = ["Half maximum (mV)","Slope","Minimum Tau (s)"];
ct=cell(1,length(ylab));
pa = zeros(1,length(ylab));
for a=1:length(axs43)
    axs43(a).YLabel.String = ylab{a};
    disp(ylab{a})
    [ct{a},pa(a)] = quickanova(Yiah(:,a),Yigh);
    makesig(ct{a},pa(a),axs43(a));
end


%atfun = @(p,x) (p(1) + p(2))./(1 
axs22 = axes('Position',[0.55 0.55 0.4 0.4]);
for n=1:length(neuron)
    Ym = mean(DKp{n}(:,:,4));
    err = nanstd(DKp{n}(:,:,4))/sqrt(size(DKp{n}(:,:,4),1));
    Ye = [Ym + err;Ym - err];
    plot(command{n}(1:2,:)+n-2,Ye,'Color',colorss(n,:));hold on
    scatter(command{n}(1,:)+n-2,Ym,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
end
% legend(names)
axs22.Title.String = 'Inactivation time constant';
axs22.YLim(1) = 0;
axs22.YLabel.String = 'Time constant (s)';
axs22.XLim(1) = -60;


subplot(axs22)
funit = @(p,x) (p(1) - p(4))./(1+exp((x - p(2))./p(3))) + p(1);
% parain = [0.25, -17,10; ...
%           0.25, -50,10; ...
%           0.25, -17,8];

%%
parain = [0.26, -14,7; ...
          0.18, -60,15; ...
          0.15, -14,8];
      
if exist('pltf','var')
    delete(pltf)
end
pltf = gobjects(1,3);
for n=1:3
    pltf(n) = plot(-60:30, funit([parain(n,:), prec(n)],-60:30),'Color',colorss(n,:)) ;
end

% pltf.YData = funit([0.25, -17,10],-60:30);
% atplt(1).YData = tfun([-13, 5, 0.06],-40:30);
%%
save('K','data','command')
%%
figure
bcolors = makecolor(0.5);
for n=1:length(neuron)
    fW2 = 4000:14000;
    X = linspace(0,length(fW2)/1e4,length(fW2));
    Y = squeeze(data{n}(fW2,1,end,1,:));
    Y = Y./repmat(max(Y(420:1e4,:)),length(fW2),1);
    
    fx = [X fliplr(X)];
    err = std(Y,[],2)/sqrt(size(Y,2));
    err = [-err';err'] + mean(Y');
    fy = [err(1,:) fliplr(err(2,:))];
    fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(n,:));hold on

    %plot(X,Y,'color',bcolors(p,:));hold on
    plot(X,mean(Y,2),'color',dcolors(n,:));hold on
end
ax = gca;
ax.YLim = [-0.05 1.05];
ax.XLim = [0 0.3];
%%
DK = squeeze(data{1}(fW2,1,:,3,:));
figure
for p=1:size(DK,2)
    plot(mean(DK(:,p,:),3)/67.03,'Color',colors(p,:));hold on% 67.03
    plot(mean(Ndata(:,p,:),3)/36.95,'Color',colors(p,:),'LineWidth',3)
end

