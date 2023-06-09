
% 
% folder = 'D:\Data\2019\';
% 
% neuron = ["B51","B64","B8"];
% info=cell(1,3);
% data = cell(1,3);
% datai = cell(1,3);
% datac = cell(1,3);
% datar = cell(1,3);
% vrange = [-60 , -20];
% command = cell(1,3);
% commandi = cell(1,3);
% xidx = [1160 4300];
% for n=1:length(neuron)
%     [~,~,info{n}] = xlsread([folder 'NaPttxAxonPP.xlsx'],neuron{n});
%     info{n} = info{n}(3:end,2:end);
%     info{n} = cellfun(@num2str,info{n},'UniformOutput',false);
% 
%     [data{n},si] = readsdata(folder,info{n}(:,[4 6]),1:2e4,4.5e3,10, 2:4, true);
%     
%     
%     command{n} = round(squeeze(data{n}(5000,2,:,1,:))'/10)*10;
%     if any(command{n}(:) ~= repelem(command{n}(1,:)',size(command{n},1))); error('pulses different'); end
%     sel = command{n}(1,:)>=vrange(1) & command{n}(1,:)<=vrange(2);
%     command{n} = command{n}(:,sel);
%     data{n}(:,:,~sel,:,:) = [];
%     data{n} = data{n}(:,:,:,3,:);
%     
%     dta = readpass(folder,info{n}(:,1:3),false);
%     datac{n} = dta(:,:,:,1);
% end
% 
% 
% 
% save('Napp.mat','data','datac','command','info')

disp('>>>> finished reading data <<<<<')

%% load data

neuron = ["B51","B64","B8"];

load('Napp.mat','data','datac','command','info');
%% passive properties
disp('passive')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
pass = cell(1,3);
cfun = @(p,x) p(1).*exp(-x./p(2)) + p(3);
cp0 = [20,0.2,-80];
X = linspace(0,1,5000);
for n=1:length(neuron)
    for e=1:size(info{n},1)
        cparam = lsqcurvefit(cfun,cp0,X',datac{n}(1:5000,1,e,1),...
                    [-inf,0,-inf],[inf,inf,inf],opts);
        pass{n}(e,1:3) = cparam;

        pass{n}(e,4) = pass{n}(e,1)/range(datac{n}(:,2,e,1));
        pass{n}(e,5) = pass{n}(e,2)/pass{n}(e,4);
        pass{n}(e,6) = pass{n}(e,5)*100;
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
        ax = subplot(3,4,e);
        for t=1:size(data{n},3)
            plot(squeeze(data{n}(:,1,t,s,e)),'Color',colors(t,:));hold on
        end
        ax.YLim = [-30, 20];
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

% fits
disp('run fits')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
fW = 3900:10300;%4090:14300
X = linspace(0,length(fW)/1e4,length(fW));
DKp = cell(1,3);
ofst = 0;%60;
for n=1:length(neuron)
    measin{n} = squeeze(mean(data{n}(14500:14520,1,:,:,:)));% 14200:14300
    init = mean(data{n}(14500:14520,1,1,:,:));
    init = repmat(init,[1,1,size(data{n},3),1,1]);
    measin{n} = ones(size(measin{n})) - measin{n}./squeeze(init);%measin{n}(t,s,e)
    measin{n}(measin{n}<0) = 0;
    for e=1:size(info{n},1)
        disp([neuron{n},'  ',info{n}{e,1}])
        for t=1:size(data{n},3)
            yy = data{n}(fW(1+ofst:end),1,t,s,e);
            fun = @(p,x) p(1).*(1 - exp(-x./0.001)).*(1 - 0.40.*(1-exp(-x./p(4))));
            p0 = [ mean(yy(100:500))   , 0.001, 0.40,  0.5];
            lb = [-30 -inf 0 0.2]; 
            ub = [0 inf 1 5];
            DKp{n}(e,t,:,s) = lsqcurvefit(fun,p0 ,X(1+ofst:end)',yy,lb,ub,opts);
        end
    end
end
%%

fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
KAa = cell(1,3);
maxg = cell(1,3);
for n=1:3
    gg = squeeze(DKp{n}(:,:,1,:)) ./repmat((command{n} - 60),[1,1]);
    gg(gg>10) = 0;
    m = max(gg,[] ,2);
    maxg{n} = m;
    m = repmat(m,1,size(command{n},2));
    KAa{n} = gg./m;%repmat(gg(:,end,:),1,size(gg,2),1);%m;%KAp{n}(:,:,1)./(m.*(command{n} + 70));
end

% check fits

X = linspace(0,length(fW)/1e4,length(fW));
fun = @(p,x) p(1).*(1 - exp(-x./p(2))).^2.*(1 - p(3).*(1 - exp(-x./p(4))));
for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} 'fits']))% 
    figure('Name',[neuron{n} 'fits'],'NumberTitle','off','Position',[(n-1)*600 100 600 900]);
    for e=1:size(info{n},1)
        axp = subplot(3,4,e);
        plt = gobjects(1,size(data{n},3));
        for t=1:size(data{n},3)
            plt(t) = plot(X,data{n}(fW,1,t,s,e),'Color',colors(t,:));hold on
            plot(X,fun(DKp{n}(e,t,:,s),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-16 5];
    end
    legend(plt,num2str(command{n}(1,:)'))
end


%% scatters
disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
fig1 = figure('Name','scatter','NumberTitle','off','Position',[00 100 1200 900],'Color','w');
names = ["B51","B64","B8"];

% activation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
axs0 = axes('Position',[0.05 0.55 0.3 0.4]);
colors = makecolor;

bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2))).^1;
bp0 = [-50, 1];
KAap = cell(1,3);
indiv = 0;

for n=1:length(neuron)
    plti = gobjects(size(KAap{n},1));

    Yy = KAa{n}(:,:,s);
    iqr = prctile(Yy,[25 75]);
    err = nanstd(Yy)/sqrt(size(KAa{n},1));
    err = [-err;err] + nanmean(Yy);

    Xx = command{n}(1,:)+n-2;
    for e=1:size(DKp{n},1)
        isn = ~isnan(KAa{n}(e,:,s));
        Xc = command{n}(e,isn);
        Y = Yy(e,isn);
        Xcf = Xc(~isinf(Y));
        Yf = Y(~isinf(Y));
        KAap{n}(e,:) = lsqcurvefit(bfun,bp0,Xcf,Yf,[-inf,-inf],[inf,inf],opts);
        if indiv
            plti(e) = plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
            scatter(command{n}(e,:)+n,KAa{n}(e,:,3),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
        end
    end
    cclr = 'k';
    if ~indiv
        cclr = colorss(n,:);
    end
    scatter(Xx,nanmean(Yy),20,'d','MarkerFaceColor',cclr,'MarkerEdgeColor',cclr);hold on
    plot(repmat(Xx,2,1),err,'Color',cclr);hold on
    aplt(n) = plot(-80:20,bfun(mean(KAap{n}),-80:20),'Color',cclr);hold on      
end
axs0.Title.String = 'Activation Curve';
axs0.YLabel.String = 'A';
axs0.XLabel.String = 'Command potential (mV)';
axs0.XLim = [-80,20];
if indiv
    legend(plti,string(info{n}(:,1)),'Location','southeast')
end

%
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
    Ys = DKp{n}(:,end,1,s)/(command{n}(1,end) - 60);%./pass{n}(:,6);
    YaA = [YaA; Ys];
    scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Ys,[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));
end 
axs1(1).YLabel.String = 'Half activation (mV)';
axs1(2).YLabel.String = 'Slope';
% axs1(3).YLabel.String = 'Current Density (nA/mm^2)';
axs1(3).YLabel.String = 'Maximum Conductance (g)';
axs1(3).YLim(1) = 0;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

disp('Half activation')
quickanova(Yah(:,1),Ygh);
disp('Slope')
quickanova(Yah(:,2),Ygh);
disp('Current Density')
quickanova(YaA,Ygh);


%

% inactivation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
axs2(1) = axes('Position',[0.05 0.05 0.3 0.4]);
% axes(axs0)

bfun = @(p,x) (1 - p(3))./(1+exp((p(1) - x)/p(2))) + p(3);
bp0 = [-60, -1,0];
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
%     iaplt(n) = plot(-80:20,bfun(mean(KAiap{n}),-80:20),'Color',colorss(n,:),'LineStyle','--');hold on
    iaplt(n) = plot(-80:20,bfun([-62, -1, 0.6],-80:20),'Color',colorss(n,:),'LineStyle','--');hold on
end
legend(aplt,names,'Location','southeast')
legend(iaplt,names,'Location','southeast')
axs2(1).Title.String = 'Inactivation Curve';
axs2(1).YLabel.String = 'B';
axs2(1).XLabel.String = 'Command potential (mV)';
axs2(1).YLim = [0 1];
axs2(1).Box = 'off';
axs2(1).XLim = [-80,-9];


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
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);
disp('Half activation')
quickanova(Yiah(:,1),Yigh);
disp('Slope')
quickanova(Yiah(:,2),Yigh);
disp('Maximum Inactivation')
quickanova(Yiah(:,3),Yigh);

%%

close(findobj(0,'Name', 'scatter time constant'))% 
figure('Name','scatter time constant','NumberTitle','off');




indv = 0;

%atfun = @(p,x) (p(1) + p(2))./(1 
axs22 = axes('Position',[0.35 0.35 0.6 0.6]);
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
funit = @(p,x) (p(1) - p(4))./(1+exp((x - p(2))./p(3))) + p(4);
p0 = [1 -45 5 0.001];
itp = cell(1,3);
for n=1:length(neuron)
    
    Ym = mean(DKp{n}(:,:,4,s));
    err = nanstd(DKp{n}(:,:,4,s))/sqrt(size(DKp{n}(:,:,4,s),1));
    Ye = [Ym + err;Ym - err];
    cclr = colorss(n,:);
    if ~indiv
        plot(command{n}(1:2,:)+n-2,Ye,'Color',cclr);hold on
        scatter(command{n}(1,:)+n-2,Ym,20,'d','MarkerFaceColor',cclr,'MarkerEdgeColor',colorss(n,:));hold on
    end
    if indv
        close(findobj(0,'Name', ['individual time constant ' neuron{n}]))% 
        figure('Name',['individual time constant ' neuron{n}],'NumberTitle','off');
        for e=1:size(DKp{n},1)
            subplot(4,3,e)
            Yi = DKp{n}(e,:,4,3);
            xi = command{n}(1,:)+n-2;
            scatter(xi,Yi,20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
            itp{n}(e,:) = lsqcurvefit(funit,p0,xi,Yi,[-inf,-inf,-inf,-inf],[inf,inf,inf,inf],opts);
            xfi = min(xi):max(xi);
            plot(xfi,funit(itp{n}(e,:),xfi),'Color',colors(e,:));hold on
            ax = gca;
            ax.YLim = [0 5];
        end
    end
end


xe = -65:0.5:0;
% funit = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))./p(4))) + p(2);
% plot(xe,funit([2 0.25 -45 2],xe),'k:');hold on

funit = @(p,x) (p(1) - p(2))./((1+exp((x - p(3))./p(4))).*(1+exp((x - p(5))./p(6)))) + p(2);
plot(xe,funit([2.5 0.25 -55 -1.5 -42 1.5],xe),'k');hold on


ax = gca;
ax.YLim = [0 3];

% % legend(names)
% axs22.Title.String = 'Inactivation time constant';
% axs22.YLim(1) = 0;
% axs22.YLabel.String = 'Time constant (s)';
% axs22.XLim(1) = -60;

%%

% %%

