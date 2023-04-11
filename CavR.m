% folder = '\\uthsch-nas.uthouston.edu\ms_nba\ByrneLab\cneveu\Data\2019\';
if strcmp(getenv('COMPUTERNAME'),'MSNBA160674') || strcmp(getenv('COMPUTERNAME'),'MSNBA-C000045')
%     folder = 'C:\Users\cneveu\Desktop\Data\2019\';
%     folder = 'F:\Data\2019\';
    folder = 'C:\Users\cneveu\My Drive\Data\2019\';
end


neuron = ["B51","B64","B8"];
info=cell(1,3);
data = cell(1,3);
command = cell(1,3);
commandi = cell(1,3);
datac = cell(1,3);
datai = cell(1,3);
mv = -40:10:20;
for n=1:length(neuron)
%     [~,~,info{n}] = xlsread([folder 'Ltype.xlsx'],[neuron{n},'_old']);
    [~,~,info{n}] = xlsread([folder 'Rtype.xlsx'],neuron{n});
    info{n} = info{n}(2:end,2:end);
    info{n} = cellfun(@num2str,info{n},'UniformOutput',false);

    [data{n},si,command{n}] = readsdata(folder,info{n}(:,1:2),1:1e4,4.5e3);

    data{n}(:,[1 3],:,3,:) = data{n}(:,[1 3],:,3,:) - repmat(data{n}(:,[1 3],1,3,:),1,1,size(data{n},3),1,1);% subtracts -40 current trace
    data{n}(:,[1 3],:,3,:) = lowpassf(data{n}(:,[1 3],:,3,:),'Fpass',300,'Fstop',500,'Fs',1e6/si);
    
    
    mvidx = arrayfun(@(x) find(nanmean(command{n},1)==x),mv);
    data{n} = data{n}(:,:,mvidx,:,:);
    command{n} = command{n}(:,mvidx);
    
    if n>1
        dta = readpass(folder,info{n}(:,4:6)); 
    else
        dta = readpass(folder,info{n}(:,3:5));
    end
    datac{n} = dta(:,:,:,1);
    
    % inactivation
    logi = ~ismember(info{n}(:,3),{'NaN'});
    if any(logi)
        if n>1
            [id,si] = readsdata(folder,info{n}(logi,3),1.6e4:1.999e4,1.6e4,10);%1.6e4
            datai{n}(:,:,:,:,logi) = id;
            datai{n}(:,[1 3],:,3,logi) = lowpassf(id(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
            commandi{n} = round(squeeze(datai{n}(find(logi,1),2,:,1,:))'/10)*10;
        else
            [~,~,info51] = xlsread([folder 'Ltype.xlsx'],'B51');
            info51 = info51(2:end,2:end);
            info51 = cellfun(@num2str,info51,'UniformOutput',false);
            [id,si] = readsdata(folder,info51(:,7),1.6e4:1.999e4,1.6e4,10);%1.6e4
            datai{n} = id;
            datai{n}(:,[1 3],:,3,:) = lowpassf(id(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
            commandi{n} = round(squeeze(datai{n}(1,2,:,1,:))'/10)*10;
        end
    end
%     cont = mean(datai{n}(2200:2300,[1 3],end,:,:));
%     cont = repmat(cont,1,1,size(datai{n},3),1,1);
    cont = repmat(datai{n}(:,[1 3],end,:,:),1,1,size(datai{n},3),1,1);
    datai{n}(:,[1 3],:,:,:) = datai{n}(:,[1 3],:,:,:) - cont;
end

xpos = [0,600,1200];
ypos = [0,250,500];
small = [600,420];
big = [1000,9000];

disp('>>>> finished reading data <<<<<')

%data{1}(:,[1 3],end,:,3) = nan;

%% Activation current traces
s=3;
colors = makecolor;
dcolors = makecolor(-0.3);
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtraction']))% 
    fig = figure('Name',[neuron{n} ' subtraction'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(3,4,e);
        for t=1:size(data{n},3)
            plot(squeeze(data{n}(:,2,t,3,e)),'Color',colors(t,:));hold on
        end
        ax.YLim = [-20, 5];
        ax.Title.String = info{n}{e,1};
        if e==1
            legend(num2str(command{n}(1,:)'),'Location','SouthWest')           
        end
    end
end

%% Inactivation current traces
s=3;
colors = makecolor;
dcolors = makecolor(-0.3);
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' inac subtr']))% 
    fig = figure('Name',[neuron{n} ' inac subtr'],'NumberTitle','off');
    for e=1:size(datai{n},5)
        ax = subplot(3,4,e);
        for t=1:size(datai{n},3)
            plot(squeeze(datai{n}(:,1,t,s,e)),'Color',colors(t,:));hold on
        end
        ax.YLim = [-50, 10];
        ax.Title.String = info{n}{e,1};
        if e==1
            legend(num2str(commandi{n}(1,:)'),'Location','SouthWest')           
        end
%         pline(ax,[2200 2300]);
    end
end

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
% 
%         ax = subplot(2,3,e);
%         p2 = plot(X,squeeze(datac{n}(:,2,e,1)));hold on
%         pp = plot(X,cfun(cparam,X),'k');hold on
%         uistack(p2,'bottom')
%         ax.YLim = [-140, -70];
%         ax.Title.String = info{n}{e,7};
    end
end

%% fit inactivation
disp('fit inact')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
ifun = @(p,x) p(1).*exp(-x./p(2)) + p(3);
p0 = [0,500,0];
Cavi = cell(1,3);
Cavif = cell(1,3);
for n = 1:length(neuron)
    for e = 1:size(datai{n},5)
        for t = 1:size(datai{n},3)
            y = squeeze(datai{n}(700:2200,3,t,3,e));
            fit = lsqcurvefit(ifun,p0,1:length(y),y',[-50,100,-Inf],[0,Inf,Inf],opts);
            Cavi{n}(t,e,:) = fit;
            Cavif{n}(t,e) = ifun(fit,-350);
        end
%         Cavif{n}(:,e) = Cavif{n}(:,e)-Cavif{n}(end,e);
    end
end


for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} 'fits inact']))% 
    figure('Name',[neuron{n} 'fits inact'],'NumberTitle','off');
    for e = 1:size(datai{n},5)
        axp = subplot(3,4,e);
        plt = gobjects(1,size(datai{n},3));
        for t=1:size(datai{n},3)
            Y = datai{n}(350:2200,1,t,s,e);
            X = (-350:2200-700)';
            plt(t) = plot(X,Y,'Color',colors(t,:));hold on
            plot(X,ifun(Cavi{n}(t,e,:),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-50 10];
        axp.XLim = [-350 1500];
    end
    legend(plt,num2str(commandi{n}(1,:)'))
end
%%
tifun = @(p,x) (0 - p(3))./(1+exp((p(1) - x)/p(2))) + p(3);
p0 = [-30,5,-10];
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];
SSin = cell(1,3);
SSinMax = cell(1,3);  % very sloppy to adjust for SSin for boxplots.
for n=1:length(neuron)
    close(findobj(0,'Name', ['inascatter' neuron{n}]))
    fig = figure('Name',['inascatter' neuron{n}],'NumberTitle','off');
    cstp = (1 - max(colorss(n,:)))/size(datai{n},5);
    xx = commandi{n}(1,:);
    cnt = 0;
    for e=1:size(datai{n},5)
        yy = Cavif{n}(:,e)';
        if n==3 && e==5
            yy(end) = 0;
        end
%         SSin{n}(e,:) = lsqcurvefit(tifun,p0,[xx,ones(1,20)*10],[yy,zeros(1,20)],[-Inf,0,-Inf],[0,Inf,0],opts);
        SSin{n}(e,:) = lsqcurvefit(tifun,p0,xx,yy,[-Inf,0,-Inf],[0,Inf,0],opts);
        SSinMax{n}(e) = SSin{n}(e,3);
        yy = -yy/SSin{n}(e,end) + 1;
        SSin{n}(e,2:3) = [ -SSin{n}(e,2) , 1];% transpose curve
        subplot(2,3,e) 
        scatter(xx(:),yy(:),'d','MarkerFaceColor',colorss(n,:)+cnt*cstp,'MarkerEdgeColor',colorss(n,:)+cnt*cstp);hold on
%         plot(xx(:),yy(:),'Color',colorss(n,:)+e*0.15);hold on
        plot(-70:10,tifun([SSin{n}(e,1:2) 1],-70:10),'Color',colorss(n,:)+cnt*cstp);hold on
        cnt = cnt+1;
    end
end


%% fits
disp('run fits')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
% measin = cell(1,3);
% fW = 4320:9320;
fW = 4390:5590;%4340:5540%4320:6720;
DKp = cell(1,3);
crop = [ 80 ,20, 80];%[40,1,40]
fixt = [0.010, 0.02, 0.07];%[0.03 , 0.02, 0.07]
for n=1:length(neuron)
    for e=1:size(info{n},1)
        disp([neuron{n},'  ',info{n}{e,1}])
        for t=1:size(data{n},3)
            X = linspace(0,length(fW)/1e4,length(fW));
            yy = data{n}(fW,3,t,s,e);
%             if command{n}(e,t)==20 && n==1 && e==4
%                 yy(1737:end) = [];%yy(1737:3575) = [];
%                 X(1737:end) = [];%X(1737:3575) = [];
%             elseif command{n}(e,t)==20 && n==2 && e==1
%                 yy(1871:end) = [];
%                 X(1871:end) = [];
%             end
            
            if isnan(yy(1))  
                DKp{n}(e,t,:) = nan;continue;
            end
%             if n>1
            inn = tifun(mean(SSin{n}),command{n}(e,t));
            fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - inn.*(1-exp(-x./p(4))));
            
            p0 = [ min(yy(crop(n):end)) , 0.01  , inn, 0.02];
%             fun = @(p,x) p(1).*(1 - exp(-x./fixt(n))).*(1 - inn.*(1-exp(-x./p(4))));
%             p0 = [ min(yy(crop(n):end)) , fixt(n) , inn, 0.1];
%             else
%                 fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - p(3).*(1-exp(-x./p(4))));
%                 p0 = [ max(yy(crop:end)) , 0.01  , 0, 0.5];
%             end
%                 DKp{n}(e,t,:,s) = lsqcurvefit(fun,p0,X(75:end)',yy,[max(yy),0,0,0.25],[max(yy)*1.5,0.25,1,Inf],opts);
            DKp{n}(e,t,:) = lsqcurvefit(fun,p0,X(crop(n):end)',yy(crop(n):end),[-100,0.001,0,0.025],[0,0.05,1,0.5],opts);
        end
    end
end
% mean(DKp{1}(:,end,1)/(command{1}(1,end)-60))) to get max conductance

fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - p(3).*(1 - exp(-x./p(4))));
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} 'fits']))% 
    figure('Name',[neuron{n} 'fits'],'NumberTitle','off','Position',[xpos(n), ypos(1), 600,420]);
    for e=1:size(info{n},1)
        axp = subplot(2,3,e);
        plt = gobjects(1,size(data{n},3));
        for t=1:size(data{n},3)
            plt(t) = plot(X,data{n}(fW,1,t,s,e),'Color',colors(t,:));hold on
            plot(X,fun(DKp{n}(e,t,:),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-20 10];
        pline(axp,crop(n)/1e4);
    end
    legend(plt,num2str(command{n}(1,:)'))
end

fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
KAa = cell(1,3);
maxg = cell(1,3);
for n=1:3
    gg = squeeze(DKp{n}(:,:,1,:)) ./ (command{n} -60);%repmat((command{n} -60),[1,1]);
    gg(gg>10) = 0;
    m = max(gg,[] ,2);
    maxg{n} = m;
    m = repmat(m,1,size(command{n},2));
    KAa{n} = gg./m;%KAp{n}(:,:,1)./(m.*(command{n} + 70));
end






%% scatters

disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
figure('Name','scatter','NumberTitle','off','Position',[1000 100 1000 800]);


% activation
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
axs0 = axes('Position',[0.05 0.55 0.3 0.4]);
colors = makecolor;
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];

bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2)));
bp0 = [-10, 1];
KAap = cell(1,3);
aplt = gobjects(1,length(neuron));
for n=1:length(neuron)
%     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    Yy = KAa{n}(:,:);
%     Yy = Yy - repmat(min(KAa{n},[],2),1,size(Yy,2));
    iqr = prctile(Yy,[25 75]);
    err = nanstd(Yy,0,1)/sqrt(size(KAa{n},1));
    err = [-err;err] + nanmean(Yy,1);

    Xx = command{n}(1,:)+n-2;
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    scatter(Xx,nanmean(Yy,1),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    for e=1:size(DKp{n},1)
        isn = ~isnan(KAa{n}(e,:));
        Xc = command{n}(e,isn);
        Y = Yy(e,isn);
        KAap{n}(e,:) = lsqcurvefit(bfun,bp0,Xc,Y,[-50,0],[10,20],opts);
%         plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%         scatter(command{n}(e,:)+n,KAa{n}(e,:),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    aplt(n) = plot(-80:20,bfun(nanmean(KAap{n},1),-80:20),'Color',colorss(n,:));hold on
end
legend(aplt,neuron,'Location','southeast')
axs0.Title.String = 'Activation Curve';
axs0.YLabel.String = 'Activation';
axs0.XLabel.String = 'Command potential (mV)';
axs0.XLim = [-80,21];
axs0.TickDir = 'out';


%boxplots
axs1(1) = axes('Position',[0.4 0.55  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.55 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.55  0.12 0.4]);
Yah = zeros(0,2);
Ygh = cell(0,1);
YaA = zeros(0,2);%DKp{n}(e,t,:,s)
for n=1:length(neuron)
    Yah = [Yah;KAap{n}]; %#ok<AGROW>
    Ygh = [Ygh;repmat({neuron{n}},size(KAap{n},1),1)];
    for p=1:2
        axes(axs1(p))
        scatter(ones(size(KAap{n},1),1)*n,KAap{n}(:,p),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(KAap{n}(:,p),[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
    axes(axs1(3))
    Ys = DKp{n}(:,end,1)/(command{n}(1,end) - 60);%./pass{n}(:,6);
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
set(axs1,'XTickLabel',neuron);
set(axs1,'XLim',[0.2 3.8]);

disp('Half activation')
[cta{1},paa(1)] = quickanova(Yah(:,1),Ygh);
makesig(cta{1},paa(1),axs1(1));
disp('Slope')
[cta{1},paa(1)] = quickanova(Yah(:,2),Ygh);
makesig(cta{1},paa(1),axs1(2));
disp('Current Density')
[cta{1},paa(1)] = quickanova(YaA,Ygh);
makesig(cta{1},paa(1),axs1(3));




% inactivation
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
% axs2(1) = axes('Position',[0.05 0.05 0.3 0.4]);
axes(axs0)

bfun = @(p,x) (1 - p(3))./(1+exp((p(1) - x)/p(2))) + p(3);
bp0 = [-10, -10,0];
KAiap = cell(1,3);%DKp{n}(e,t,:,s)
iaplt = gobjects(1,3);
for n=1:length(neuron)
    Yi = -Cavif{n}./SSinMax{n}+1;
    Yi = Yi';
%     if n==1
%         Yi = squeeze(DKp{n}(:,3:end,3));
%         cmv = command{n}(:,3:end); 
%     end
%     scatter(command{n}(:)+n-2, 1-Yi(:), 20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    iqr = prctile(Yi,[25 75]);
    err = nanstd(Yi,0,1)/sqrt(size(Yi,1));
    err = [-err;err] + 1-nanmean(Yi,1);

    Xx = commandi{n}(1,:)+n-2;

    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    scatter(Xx,1-nanmean(Yi,1),20,'d','MarkerFaceColor','w','MarkerEdgeColor',colorss(n,:));hold on

    SSin{n}(:,2) = abs(SSin{n}(:,2));
    iaplt(n) = plot(-80:20,tifun(mean(SSin{n}),-80:20),'Color',colorss(n,:),'LineStyle','--');hold on
end
%legend(iaplt,names,'Location','SouthEast')
% axs2(1).Title.String = 'Inactivation Curve';
% axs2(1).YLabel.String = 'B';
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
    Yiah = [Yiah;SSin{n}]; %#ok<AGROW>
    Yigh = [Yigh;repmat({neuron{n}},size(SSin{n},1),1)];
    for p=1:3
        axes(axs1(p))
        Y = SSin{n}(:,p);
        if p==3; Y = SSinMax{n};end
        scatter(ones(size(SSin{n},1),1)*n,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
axs1(1).YLabel.String = 'Half inactivation (mV)';
axs1(1).YLim = [-40 20];
axs1(2).YLabel.String = 'Slope';
axs1(3).YLabel.String = 'Maximum Inactivation';
axs1(3).YLim = [0 1];
set(axs1,'TickDir','out')
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',neuron);
set(axs1,'XLim',[0.2 3.8]);

disp('Half inactivation')
[cti{1},pai(1)] = quickanova(Yiah(:,1),Yigh);
lns = makesig(cti{1},pai(1),axs1(1));
disp('Inact Slope')
[cti{2},pai(2)] = quickanova(Yiah(:,2),Yigh);
makesig(cti{2},pai(2),axs1(2));

set(findobj('Type','axis'),'TickDir','out')
%%
% fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - p(3).*(1-exp(-x./p(4))));
% close(findobj(0,'Name', 'curve'))% 
% figure('Name','curve','NumberTitle','off');
% bcolors = makecolor(0.5);
% mvm = 20;
% cplt = gobjects(1,3);
% for n=1:length(neuron)
%     fW2 = 4452:6742;
%     X = linspace(0,length(fW2)/1e4,length(fW2));
% %     Y = squeeze(data{n}(fW2,3, end ,3,:))./repmat(pass{n}(:,6)',length(fW2),1);
%     idx = command{n}(1,:)==mvm;
%     Y = squeeze(data{n}(fW2,3,idx ,3,:));
%     Y = -Y./repmat(min(Y),length(fW2),1);
% %     Y = Y./repmat(max(Y(420:1e4,:)),length(fW2),1);
%     
%     fx = [X fliplr(X)];
%     err = nanstd(Y,[],2)/sqrt(size(Y,2));
%     err = [-err';err'] + nanmean(Y');
%     fy = [err(1,:) fliplr(err(2,:))];
%     fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(n,:));hold on
% 
%     %plot(X,Y,'color',bcolors(p,:));hold on
%     cplt(n) = plot(X,mean(Y,2),'color',dcolors(n,:));hold on
% %     Yf = fun(mean(DKp{n}(:,idx,:)),X+(fW2(1)-fW(1))*1e-4);
% %     Yf = -Yf./min(Yf);
%     Yf = zeros(length(fW2),size(info{n},1));
%     for e=1:size(info{n},1)
%         Yf(:,e) = fun(DKp{n}(e,idx,:),X+(fW2(1)-fW(1))*1e-4);
%     end
%     Yf = -Yf./repmat(min(Yf),length(fW2),1);
%     plot(X,mean(Yf,2),'LineStyle','--','color',dcolors(n,:));hold on
% end
% ax = gca;
% % ax.YLim = [-0.05 1.05];
% ax.XLim = [0 0.2];
% ax.YLim = [-1.1 0.3];
% ax.YLabel.String = 'Current (nA)';
% ax.XLabel.String = 'Time (s)';
% lg = legend(cplt,neuron,'Location','SouthEast');
% text(lg.Position(1)+0.05,lg.Position(2)+0.1,['step to ',num2str(mvm),' mV'],'Units','normalized','HorizontalAlignment','left')

%%
fun = @(p,x) p(1).*(1 - exp(-x./p(2))).*(1 - p(3).*(1-exp(-x./p(4))));
close(findobj(0,'Name', 'curve steps'))% 
figure('Name','curve steps','NumberTitle','off');
bcolors = makecolor(0.5);
cplt = gobjects(1,3);
fW2 =  4390:5590;%4340:5540%4452:6742;
X = linspace(0,length(fW2)/1e4,length(fW2));
mvm = 0:10:20;
for n=1:length(neuron)
    ax = subplot(1,3,n);
    ftline = gobjects(1,length(mvm));
    for t=1:length(mvm)
    %     Y = squeeze(data{n}(fW2,3, end ,3,:))./repmat(pass{n}(:,6)',length(fW2),1);
        idx = command{n}(1,:)==mvm(t);
        Y = squeeze(data{n}(fW2,3,idx ,3,:));
        Y = -Y./repmat(min(Y(100:end,:)),length(fW2),1);
        Y = -Y/min(mean(Y(100:end,:),2));
    %     Y = Y./repmat(max(Y(420:1e4,:)),length(fW2),1);

        fx = [X fliplr(X)];
        iqr = [prctile(Y',25),fliplr(prctile(Y',75))];
        err = nanstd(Y,[],2)/sqrt(size(Y,2));
        err = [-err';err'] + repmat(nanmean(Y'),2,1);
        fy = [err(1,:) fliplr(err(2,:))];
        fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(t,:));hold on

        %plot(X,Y,'color',bcolors(p,:));hold on
        cplt(t) = plot(X,mean(Y,2),'color',dcolors(t,:));hold on
    %     Yf = fun(mean(DKp{n}(:,idx,:)),X+(fW2(1)-fW(1))*1e-4);
    %     Yf = -Yf./min(Yf);
        Yf = zeros(length(fW2),size(info{n},1));
        for e=1:size(info{n},1)
            Yf(:,e) = fun(DKp{n}(e,idx,:),X);%+(fW2(1)-fW(1))*1e-4);
        end
        Yf = -Yf./repmat(min(Yf),length(fW2),1);
        ftline(t) = plot(X,mean(Yf,2),'LineStyle','--','color',dcolors(t,:),'LineWidth',2);hold on
    end
    uistack(ftline,'top')
    ax.XLim = [0 0.12];
    ax.YLim = [-1.1 0];
    ax.YLabel.String = 'Current (nA)';
    ax.XLabel.String = 'Time (s)';
    lg = legend(cplt,num2str(mvm'),'Location','NorthWest');
end


%%

close(findobj(0,'Name', 'scatter time constant'))% 
figure('Name','scatter time constant','NumberTitle','off','Position',[0 50 1800 850]);


atp = cell(1,3);
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
mint = 0.001;
for n=1:length(neuron)
    for e=1:size(DKp{n},1)
        Yd = DKp{n}(e,2:end,2);
        xt = command{n}(e,2:end);
        tfun = @(p,x) (p(1) - mint)./(1+exp((x - p(3))/-p(4)))./(1+exp((x - (p(3) + p(5)))/p(4))) + mint;
        tp0 = [0.1, mint, -20, 8, 20];%max(Yd)
        atp{n}(e,:) = lsqcurvefit(tfun,tp0,xt,Yd,[0,0,-inf,-inf,0],[inf,inf,inf,inf,inf],opts);
        atp{n}(e,5) = sum(atp{n}(e,[3 5]));
    end
end
tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(6))) + p(2);

% activation time constant scatter
axs21 = axes('Position',[0.05 0.55 0.3 0.4]);
atplt = gobjects(1,3);
for n=1:length(neuron)
    Yd = DKp{n}(:,2:end,2);
    Ym = mean(Yd);
    err = nanstd(DKp{n}(:,2:end,2))/sqrt(size(DKp{n}(:,2:end,2),1));
    Ye = [Ym + err;Ym - err];
    xt = command{n}(:,2:end);
%     scatter(xt(:),Yd(:),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    plot(command{n}(1:2,2:end)+n-2,Ye,'Color',colorss(n,:));hold on
    scatter(command{n}(1,2:end)+n-2,Ym,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
%     atplt(n) = plot(-80:30,tfun(median(atp{n}),-80:30),'Color',colorss(n,:));hold on
    if n==1
        atplt(n) = plot(-80:30,tfun([0.03 0.001 -10 5.3 14 5.6],-80:30),'Color',colorss(n,:));hold on
    elseif n==2
        atplt(n) = plot(-80:30,tfun([0.02 0.001 13 12.7 25  13],-80:30),'Color',colorss(n,:));hold on
    else
        atplt(n) = plot(-80:30,tfun([0.076 0.001 -1 7.2 23  7.2],-80:30),'Color',colorss(n,:));hold on
    end
end
% atplt(1) = plot(-80:30,tfun([0.035    0.0005   -2.5    5.6    7.75],-80:30),'Color',colorss(1,:));hold on
% % atplt(1) = plot(-80:30,tfun([0.035    0.0005   -4    7    7.75],-80:30),'Color',colorss(1,:));hold on
% atplt(2) = plot(-80:30,tfun([0.014    0.0005   0         8    18],-80:30),'Color',colorss(2,:));hold on
% atplt(3) = plot(-80:30,tfun([0.12    0.0005   7         7    16.5],-80:30),'Color',colorss(3,:));hold on
% legend(atplt,neuron)
axs21.Title.String = 'Activation time constant';
axs21.XLim = [-80 30];
axs21.YLabel.String = 'Time constant (s)';




%boxplots
axs43(1) = axes('Position',[0.39 0.55  0.09 0.4]);
axs43(2) = axes('Position',[0.51 0.55  0.09 0.4]);
axs43(3) = axes('Position',[0.63 0.55  0.09 0.4]);
axs43(4) = axes('Position',[0.75 0.55  0.09 0.4]);
axs43(5) = axes('Position',[0.87 0.55  0.09 0.4]);
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
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
set(axs43,'XTick',1:3);
set(axs43,'XTickLabel',neuron);
set(axs43,'XLim',[0.2 3.8]);

ylab = ["Maximum Tau (s)","Minimum Tau (s)","h1 (mV)","Slope","h2 (mV)"];
ct=cell(1,length(ylab));
pa = zeros(1,length(ylab));
for a=1:length(axs43)
    axs43(a).YLabel.String = ylab{a};
    disp(ylab{a})
    [ct{a},pa(a)] = quickanova(Yiah(:,a),Yigh);
    makesig(ct{a},pa(a),axs43(a));
end


opts = optimset('Display','off','Algorithm','levenberg-marquardt');
itplt = gobjects(1,3);
axs22 = axes('Position',[0.05 0.05 0.3 0.4]);
mvmf = -30:10:20;
atpi = cell(1,3);
tpi0 = [-5, 6, 0.03, 0.055];
tifun = @(p,x) (p(4) - p(3))./(1+exp((x - p(1))/p(2))) + p(3);
for n=1:length(neuron)
    iidx = arrayfun(@(x) any(mvmf==x,'all'),command{n}(1,:));
    X = command{n}(1,iidx);
    Y = DKp{n}(:,iidx,4);
    if n==1
        X = X(:,3:end);
        Y = Y(:,3:end);
    end
    
    Ym = mean(Y);
    iqr = prctile(Y,[25 75]);
    err = nanstd(Y)/sqrt(size(Y,1));
    Ye = [Ym + err;Ym - err];
    
    plot([X;X]+n-2,Ye,'Color',colorss(n,:));hold on
    itplt(n) = scatter(X+n-2,nanmean(Y),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    for e=1:size(Y,1)
        atpi{n}(e,:) = lsqcurvefit(tifun,tpi0,X,Y(e,:),[-inf,-inf,-inf,-inf],[inf,inf,inf,inf],opts);
    end
    if n==1
        % [0 5 0.038 0.08]
        ppltt(n) = plot(-81:30,tifun([-1 5 0.036 0.085],-81:30),'Color',colorss(n,:));hold on
    elseif n==3
        ppltt(n) = plot(-81:30,tifun([1 6 0.025 0.05],-81:30),'Color',colorss(n,:));hold on
    else
%         ppltt(n) = plot(-81:30,tifun(median(atpi{n}),-81:30),'Color',colorss(n,:));hold on
        ppltt(n) = plot(-81:0.2:30,tifun([0 1.3 0.025 0.036] ,-81:0.2:30),'Color',colorss(n,:));hold on
    end
end

legend(ppltt,neuron)
axs22.Title.String = 'Inactivation time constant';
% axs22.YLim = [0 0.5];
axs22.YLabel.String = 'Time constant (s)';
axs22.XLim = [-81,30];
axs22.YLim(1) = 0;

%boxplots
axs33(1) = axes('Position',[0.40 0.05  0.10 0.4]);
axs33(2) = axes('Position',[0.55 0.05  0.10 0.4]);
axs33(3) = axes('Position',[0.70 0.05  0.10 0.4]);
axs33(4) = axes('Position',[0.85 0.05  0.10 0.4]);
Yiah = zeros(0,2);
Yigh = cell(0,1);
YiaA = zeros(0,3);%DKp{n}(e,t,:,s)
for n=1:length(neuron)
    Yiah = [Yiah;atpi{n}]; %#ok<AGROW>
    Yigh = [Yigh;repmat({neuron{n}},size(atpi{n},1),1)]; %#ok<AGROW>
    for p=1:size(atpi{n},2)
        axes(axs33(p))
        Y = atpi{n}(:,p);
        scatter(ones(size(atpi{n},1),1)*n,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
set(axs33,'XTick',1:3);
set(axs33,'XTickLabel',neuron);
set(axs33,'XLim',[0.2 3.8]);

varn = ["Half","Slope","Min","Max"];
ylab = ["Half activation (mV)","Slope","Minimum Tau (s)","Maximum Tau (s)"];
ct=cell(1,length(varn));
pa = zeros(1,length(varn));
for a=1:length(axs33)
    axs33(a).YLabel.String = ylab{a};
    disp(varn{a})
    [ct{a},pa(a)] = quickanova(Yiah(:,a),Yigh);
    makesig(ct{a},pa(a),axs33(a));
end
save('CaR','data','command')

%%
% n = 1;
% close(findobj(0,'Name', ['indv' neuron{n}]))% 
% figure('Name',['indv' neuron{n}],'NumberTitle','off');
% 
% 
% Xx = command{n}(1,1:end);
% leg = strings(0,1);
% cnt = 1;
% indv = gobjects(1);
% for e=1:size(DKp{n},1)
%     scatter(Xx,1-DKp{n}(e,1:end,3),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
%     indv(cnt) = plot(-60:20,bfun(KAiap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%     leg = [leg;string(['exp ',num2str(e)])]; %#ok<AGROW>
%     scatter(commandi{n}(e,:),Cavi{n}(:,e,1)/min(Cavi{n}(:,e,1)),'k','MarkerFaceColor','k')
%     cnt = cnt +1;
% end
% 
% legend(indv,leg,'Location','SouthWest')
% axx = gca;
% axx.Title.String = 'Inactivation Curve';
% axx.YLabel.String = 'B';
% axx.XLabel.String = 'Command potential (mV)';
% axx.YLim = [0 1];
% axx.XLim = [-60,21];

%%
figure
bcolors = makecolor(0.5);
for n=1:length(neuron)
    fW2 = 4000:10000;
    X = linspace(0,length(fW2)/1e4,length(fW2));
    Y = squeeze(data{n}(fW2,1,end,3,:));
%     Y = Y./repmat(max(Y(420:1e3,:)),length(fW2),1);
    
    fx = [X fliplr(X)];
    err = std(Y,[],2)/sqrt(size(Y,2));
    err = [-err';err'] + mean(Y');
    fy = [err(1,:) fliplr(err(2,:))];
    fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(n,:));hold on

    %plot(X,Y,'color',bcolors(p,:));hold on
    plot(X,mean(Y,2),'color',dcolors(n,:));hold on
end
ax = gca;
% ax.YLim = [-0.05 2.05];
% ax.XLim = [0 0.3];
