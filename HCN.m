
% folder = 'D:\Data\B8_data-HMC\';
% 
% neuron = ["B51","B64","B8"];
% info=cell(1,3);
% data = cell(1,3);
% 
% n = 3;
% 
% vrange = [-110 , -60];
% % vrange = [-110 , -30];
% 
% [~,~,info{n}] = xlsread([folder 'HCN.xlsx'],neuron{n});
% info{n} = info{n}(2:end,2:end);
% info{n} = cellfun(@num2str,info{n},'UniformOutput',false);
% 
% [data{n},si,command{n}] = readsdata(folder,info{n}(:,1),1:1.5e5,4e4,10,[2 1]);
% 
% vcom = nanmean(command{n});
% sel = vcom>=vrange(1) & vcom<=vrange(2);
% command{n} = command{n}(:,sel);
% data{n}(:,:,~sel,:,:) = [];
% vcom(~sel) = [];
% 
% data{n}(:,1,:,1,:) = lowpassf(data{n}(:,1,:,1,:),'Fpass',120,'Fstop',300,'Fs',1e6/si);
% data{n} = data{n} - repmat(mean(data{n}(10000:11000,:,:,:,:)),length(data{n}),1,1,1,1);
% leak = max(data{n}(25500:27000,1,:,1,:));
% 
% save('HCN.mat','data','command','info','si');

%% load data

neuron = ["B51","B64","B8"];
n = 3;
load('HCN.mat','data','command','info','si');


%% Activation current traces
s=1;
colors = makecolor;
dcolors = makecolor(-0.3);
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];

% fits
disp('run fits')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
KAa = cell(1,3);
ofst = 700;
fun = @(p,x) p(1).*(1 - exp(-x./p(2))) + p(3);
p0 = [-1 100 -1];
tailfit = cell(1,3);
KAa{n} = nan(size(data{n},[5,3]));
for e=1:size(data{n},5)
    if e==1
        fW = 65000:69000;
    else
        fW = 85000:89000;
    end
    X = linspace(0,length(fW)/1e4,length(fW));
    disp([neuron{n},'  ',info{n}{e,1}])
    for t=1:size(data{n},3)
        yy = data{n}(fW(1+ofst:end),1,t,s,e);
        if ~isnan(yy(1))
            tailfit{n}(e,t,:) = lsqcurvefit(fun,p0 ,X(1+ofst:end)',yy,[-inf, 0.5, -inf] ,[0 10 0],opts);
            KAa{n}(e,t) = fun(tailfit{n}(e,t,:),0);
        end
    end
    KAa{n}(e,:) = -(KAa{n}(e,:) - max(KAa{n}(e,:),[],'omitnan'));
    KAa{n}(e,:) = KAa{n}(e,:)./max(KAa{n}(e,:),[],'omitnan');
end

% check fits

close(findobj(0,'Name', [neuron{n} 'fits']))% 
figure('Name',[neuron{n} 'fits'],'NumberTitle','off','Position',[100 100 1700 900]);
for e=1:size(data{n},5)
    if e==1
        fW = 65000:69000;
    else
        fW = 85000:89000;
    end
    axp = subplot(3,3,e);
    plt = gobjects(1,size(data{n},3));
    for t=1:size(data{n},3)
        if ~isnan(data{n}(fW(1),1,t,s,e))
            plot(X,data{n}(fW,1,t,s,e),'Color',colors(t,:));hold on
            plot(X,fun(tailfit{n}(e,t,:),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
    end
    axp.Title.String = info{n}{e,1};
%     axp.YLim = [-30 0];
    pline(gca,X(ofst));
end


%% Fit to get time constant 
s=1;
colors = makecolor;
dcolors = makecolor(-0.3);
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];

% fits
disp('run fits')
opts = optimset('Display','off');%,'Algorithm','levenberg-marquardt');
fW = 25000:55000;%:65000;
X = linspace(0,length(fW)/1e4,length(fW));
DKp = cell(1,3);
ofst = 700;
fun = @(p,x) p(1).*(1 - exp(-x./p(2))) + p(3);
p0 = [-1 100 -1];
for e=1:size(data{n},5)
    disp([neuron{n},'  ',info{n}{e,1}])
    for t=1:size(data{n},3)
        yy = data{n}(fW(1+ofst:end),1,t,s,e);
        if ~isnan(yy(1))
            DKp{n}(e,t,:) = lsqcurvefit(fun,p0 ,X(1+ofst:end)',yy,[-inf, 0.5, -inf] ,[0 10 0],opts);
        end
    end
end
DKp{n}(DKp{n}==0) = nan;
% DKp{n}(repmat(abs(DKp{n}(:,:,1))<0.01,1,1,3)) = nan;
%

KAa{n} = DKp{n}(:,:,1);
maxg = cell(1,3);
gg = squeeze(DKp{n}(:,:,1,:)) ./repmat((command{n} +40),[1,1]);
gg(gg>10) = 0;
m = nanmax(gg,[] ,2);
maxg{n} = m;
m = repmat(m,1,size(command{n},2));
KAa{n} = gg./m;%KAp{n}(:,:,1)./(m.*(command{n} + 40));
KAa{n}(KAa{n}==-inf) = 0;
% % KAa{n}(3:end,end) = 0;

% check fits

close(findobj(0,'Name', [neuron{n} 'fits TC']))% 
figure('Name',[neuron{n} 'fits TC'],'NumberTitle','off','Position',[100 100 1700 900]);
for e=1:size(data{n},5)
    axp = subplot(3,3,e);
    plt = gobjects(1,size(data{n},3));
    for t=1:size(data{n},3)
        plt(t) = plot(X,data{n}(fW,1,t,s,e),'Color',colors(t,:));hold on
        plot(X,fun(DKp{n}(e,t,:,s),X),'Color',dcolors(t,:),'LineWidth',2);hold on
    end
    axp.Title.String = info{n}{e,1};
    axp.YLim = [-30 0];
    pline(gca,X(ofst));
end
legend(plt,num2str(nanmean(command{n})'))

%%
data{n}(25000:65030,1,:,1,1) = data{n}(25000:65030,1,:,1,1) - repmat(permute(DKp{n}(1,:,3),[1 3 2]),65030-25000+1,1,1);% < fix this with permute
data{n}(25000:85030,1,:,1,2:end) = data{n}(25000:85030,1,:,1,2:end) - repmat(permute(DKp{n}(2:end,:,3),[4 3 2 5 1]),85030-25000+1,1,1);


%% scatters
disp('run scatter')

close(findobj(0,'Name', 'scatter'))% 
figure('Name','scatter','NumberTitle','off','Position',[100 100 1700 900]);
names = ["B51","B64","B8"];

% activation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
axs1 = axes('Position',[0.05 0.55 0.3 0.4]);
colors = makecolor;

bfun = @(p,x) 1./(1+exp((x - p(1))/p(2)));
bp0 = [-80, 10];
KAap = cell(1,3);

%     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
Yy = KAa{n}(:,:,s);
%     Yy = Yy - repmat(min(KAa{n}(:,:,s),[],2),1,size(Yy,2));
iqr = prctile(Yy,[25 75]);
err = nanstd(Yy)/sqrt(size(KAa{n},1));
err = [-err;err] + nanmean(Yy);

Xx = nanmean(command{n}(:,:))+n-2;
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
scatter(Xx,abs(nanmean(Yy)),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
plot(-120:-40,bfun(mean(KAap{n}),-120:-40),'Color',colorss(n,:));hold on

% legend(aplt,names,'Location','southeast')
axs1.YLim = [0 1];
axs1.Title.String = 'Activation Curve';
axs1.YLabel.String = 'Activation';
axs1.XLabel.String = 'Command potential (mV)';



%boxplots
axs1(1) = axes('Position',[0.4 0.55  0.12 0.4]);
axs1(2) = axes('Position',[0.55 0.55 0.12 0.4]);
axs1(3) = axes('Position',[0.7 0.55  0.12 0.4]);
Yah = zeros(0,2);
Ygh = cell(0,1);

Yah = [Yah;KAap{n}]; 
Ygh = [Ygh;repmat({names{n}},size(KAap{n},1),1)];
for p=1:2
    axes(axs1(p))
    scatter(ones(size(KAap{n},1),1)*n,KAap{n}(:,p),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(KAap{n}(:,p),[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
end
axes(axs1(3))
Ys = m(:,1);%./pass{n}(:,6);
scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
Yp = prctile(Ys,[25 50 75]);
rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));

axs1(1).YLabel.String = 'Half activation (mV)';
axs1(2).YLabel.String = 'Slope';
axs1(3).YLabel.String = 'Maximum Conductance';
axs1(3).YLim(1) = 0;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

% Time contants

tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(4))) + p(2);

axs21 = axes('Position',[0.05 0.05 0.3 0.4]);
atplt = gobjects(1,3);

Yd = DKp{n}(:,:,2);
Ym = nanmean(Yd);
iqr = prctile(Yd,[25 75]);
err = nanstd(DKp{n}(:,:,2))./sqrt(sum(~isnan(Yd)));
Ye = [Ym + err;Ym - err];
xt = nanmean(command{n});
%     scatter(xt(:),Yd(:),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
plot(repmat(xt,2,1)+n-2,Ye,'Color',colorss(n,:));hold on
scatter(xt+n-2,nanmean(Yd),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on

opts = optimset('Display','off');% 'Algorithm','levenberg-marquardt');
tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - (p(3)+p(5)))/p(4))) + p(2);
atp = cell(1,3);
% tp0 = [3 0.5 -80 7 0];
tp0 = [1 0.5 -80 3 0];
for e=1:size(DKp{n},1)
    isn = ~isnan(DKp{n}(e,:,2));
    Xc = command{n}(e,isn);
    Yt = DKp{n}(e,isn,2);
    if length(Xc)>4
        atp{n}(e,:) = lsqcurvefit(tfun,tp0,Xc,Yt,[0,0,-inf,0,0],[inf,inf,inf,inf,inf],opts);
        atp{n}(e,5) = sum(atp{n}(e,[3 5]));
    end
end
% tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(4))) + p(2);
tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(4))) + p(2);
% plot(-120:-30,tfun(nanmedian(atp{n}),-120:-30),'Color',colorss(n,:));hold on
plot(-120:-30,tfun([14 0.9 -73 9 -75],-120:-30),'Color',colorss(n,:));hold on

% legend(atplt,neuron)
axs21.Title.String = 'Activation time constant';
axs21.XLim = [-120 -40];
axs21.YLim(1) = 0;
axs21.YLabel.String = 'Time constant (s)';


tfun = @(p,x) (p(1) - p(2))./(1+exp((x - p(3))/-p(4)))./(1+exp((x - p(5))/p(4))) + p(2);

%boxplots
axs43(1) = axes('Position',[0.39 0.05  0.09 0.4]);
axs43(2) = axes('Position',[0.51 0.05  0.09 0.4]);
axs43(3) = axes('Position',[0.63 0.05  0.09 0.4]);
axs43(4) = axes('Position',[0.75 0.05  0.09 0.4]);
axs43(5) = axes('Position',[0.87 0.05  0.09 0.4]);
YiaA = zeros(0,3);%DKp{n}(e,t,:,s)
% ylab = ["Maximum Tau (s)","Minimum Tau (s)","Half activation (mV)","Slope"];
ylab = ["Maximum Tau (s)","Minimum Tau (s)","h1 (mV)","Slope","h2 (mV)"];
for n=1:length(neuron)
    for p=1:size(atp{n},2)
        axes(axs43(p))
        Y = atp{n}(:,p);
        scatter(ones(size(atp{n},1),1)*n,Y,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        text(ones(size(atp{n},1),1)*n+0.1,Y,string((1:length(Y))'))
        Yp = prctile(Y,[25 50 75]);
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
        axs43(p).YLabel.String = ylab{p};
    end
end 
set(axs43,'XTick',1:3);
set(axs43,'XTickLabel',neuron);
set(axs43,'XLim',[0.2 3.8]);


%% Activation current traces



close(findobj(0,'Name', [neuron{n} ' subtraction']))% 
fig = figure('Name',[neuron{n} ' subtraction'],'NumberTitle','off');
for e=1:size(info{n},1)
    ax = subplot(3,3,e);
    for t=1:size(data{n},3)
        plot(squeeze(data{n}(:,1,t,s,e)),'Color',colors(t,:));hold on
    end
    ax.YLim = [-10, 5];
    ax.Title.String = info{n}{e,1};
    if e==1
        legend(num2str(command{n}(1,:)'),'Location','SouthWest')           
    end
end

%%
save('HCN','data','command')
