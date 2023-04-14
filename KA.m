 
% % folder = 'Data\';% the 
% folder = 'D:\Data\2019\';
% 
% neuron = ["B51","B64","B8"];
% info=cell(1,3);
% data = cell(1,3);
% datai = cell(1,3);
% datac = cell(1,3);
% datar = cell(1,3);
% vrange = [-50 , 20];
% command = cell(1,3);
% commandi = cell(1,3);
% logg = '01';
% for n=1:length(neuron)
%     [~,~,info{n}] = xlsread([folder 'Atype.xlsx'],neuron{n});
%     info{n} = info{n}(3:end,2:end);
%     info{n} = cellfun(@num2str,info{n},'UniformOutput',false);
% 
%     [data{n},si] = readsdata(folder,info{n}(:,[1 4]),1:1e4,4.5e3,10);
%     data{n}(:,[1 3],:,3,:) = lowpassf(data{n}(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
%     
%     command{n} = round(squeeze(data{n}(5000,2,:,1,:))'/10)*10;
%     if any(command{n}(:) ~= repelem(command{n}(1,:)',size(command{n},1))); error('pulses different'); end
%     sel = command{n}(1,:)>=vrange(1) & command{n}(1,:)<=vrange(2);
%     command{n} = command{n}(:,sel);
%     data{n}(:,:,~sel,:,:) = [];
%     
%     % inactivation
%     [datai{n},si] = readsdata(folder,info{n}(:,[2 5]),1.15e4:1.66e4,1.15e4,5);
%     datai{n}(:,[1 3],:,3,:) = lowpassf(datai{n}(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
%     
%     commandi{n} = round(squeeze(datai{n}(1,2,:,1,:))'/5)*5;
%     
%     %passive properties
%     for e=1:size(info{n},1)
%         [datagf,sig] = abfload([folder info{n}{e,7} '.abf'],'verbose',false);
%         datagf = datagf(round(str2double(info{n}{e,8})*1e3/sig):round(str2double(info{n}{e,9})*1e3/sig),:);
%         pulse = datagf(50:end,3) - datagf(1:end-49,3);
%         pup = logg((pulse>0.055)+1);
%         pt = strfind(pup,repelem('01',47));
%         pdn = logg((pulse<-0.055)+1);
%         pe = strfind(pdn,repelem('01',47));
%         if n==3
%             short = arrayfun(@(x) any(abs(pe-x-25000)<100),pt);
%         else
%             short = arrayfun(@(x) any(abs(pe-x-5000)<100),pt);
%         end
%         long = arrayfun(@(x) any(abs(pt-x-25000)<100),pe);       
%         pt = [pe(1) , pt(short)];
%         for t=1:2
%             if t<=length(pt)
%                 datac{n}(:,:,e,t) = datagf(pt(t)+50:pt(t)+5000,:);
%             end
%         end
%     end
% end
% 
% 
% % extracting full inactivation protocol
% 
% dataif = cell(1,3);
% for n = 1:3
%     [dataif{n},si] = readsdata(folder,info{n}(:,[2 5]),1:2e4,1.15e4,5);
%     dataif{n}(:,[1 3],:,3,:) = lowpassf(dataif{n}(:,[1 3],:,3,:),'Fpass',500,'Fstop',700,'Fs',1e6/si);
% end
% 
% disp('>>>> finished reading data <<<<<')
% save('KA','datai','data','command','commandi','datac','info','si');

%% load data

vrange = [-50 , 20];
neuron = ["B51","B64","B8"];
folder = 'Data\';

load('KA','datai','data','command','commandi','datac','info','si');



%% Activation current traces

for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtractions']))% 
    fig = figure('Name',[neuron{n} ' subtractions'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(3,3,e);
        plot(squeeze(data{n}(:,1,:,3,e)))  
        ax.YLim = [-30, 200];
        ax.Title.String = info{n}{e,1};
    end
end

%% Innactivation current traces

for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtractions inact']))% 
    fig = figure('Name',[neuron{n} ' subtractions inact'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(3,3,e);
        plot(squeeze(datai{n}(:,1,:,3,e)))  
        ax.YLim = [-30, 150];
        ax.Title.String = info{n}{e,1};
    end
    legend(string(commandi{n}(1,:)))
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
fW = 4380:6380;%9180;%6304;
X = linspace(0,range(fW)*si/1e6,length(fW));
KAp = cell(1,3);
KAa = cell(1,3);
%tau = [0.0629,0.0387,0.1361];
disp('fits')
for n=1:length(neuron)
    for e=1:size(info{n},1)
        disp([neuron{n},'  ',info{n}{e,1}])
        for t=1:size(data{n},3)
            if ~isnan(data{n}(1,1,t,3,e))
                b = mean(data{n}(9100:9190,1,t,3,e));
                %fun = @(p,x) p(1).*(1 - exp(-x/tau1(n))).*exp(-x/p(3)) + b;
                fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
                p0 = [1 0.06 b];
                KAp{n}(e,t,:) = lsqcurvefit(fun,p0,X(1:length(fW)),data{n}(fW,1,t,3,e)',...
                    [0,0,-inf],[inf,inf,inf],opts); 
            else
                KAp{n}(e,t,:) = nan;
            end
        end
    end
    KAa{n}(KAa{n}==Inf) = 0;
end
%KAp{2}(1,end,:) = nan;

fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
KAa = cell(1,3);
KAg = cell(1,3);
for n=1:3
    gg = nan(size(info{n},1),size(data{n},3));
    for e=1:size(info{n},1)
        for t=1:size(data{n},3)
            gg(e,t) = fun(KAp{n}(e,t,:),-0.0253) ./(command{n}(e,t) + 70);
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
fW = 4380:6380;%9180;%6304;
fW3 = 4380-bln:6380;%6304;
X = linspace(0,range(fW)*si/1e6,length(fW));
X2 = [fliplr(-X(1:bln)) , X];


colors = makecolor;
dcolors = makecolor(-0.3);
% fun = @(p,x) p(1).*(1 - exp(-x/p(2))).*exp(-x/p(3)) + p(4);
fun = @(p,x) p(1).*exp(-x/p(2)) + p(3);
for n=1:3
    close(findobj(0,'Name', [neuron{n} 'fits']))% 
    figure('Name',[neuron{n} 'fits'],'NumberTitle','off');
    cnt = 1;
    for e=[1  4]%expp
        if e>size(KAp{n},1);continue;end
        axp = subplot(2,ceil(length(expp)/2),cnt);
        cnt=cnt+1;
        for t=1:size(data{n},3)
            plt(t) = plot(X2,data{n}(fW3,1,t,3,e),'Color',colors(t,:));hold on
            plot(X,fun(KAp{n}(e,t,:),X),'Color',dcolors(t,:),'LineWidth',2);hold on
        end
        axp.Title.String = info{n}{e,1};
        axp.YLim = [-80 350];
    end
    legend(plt,num2str(command{n}(1,:)'))
end


%% scatters
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
%     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    iqr = prctile(KAa{n},[25 75]);
    err = nanstd(KAa{n})/sqrt(size(KAa{n},1));
    err = [-err;err] + nanmean(KAa{n});
%     fx = [command{n}(1,:) fliplr(command{n}(1,:))];
%     fy = [err(1,:) fliplr(err(2,:))];
%     fill(fx,fy,colors(n,:),'EdgeColor','none','FaceAlpha',0.5);hold on

    Xx = command{n}(1,:)+n-2;
    scatter(Xx,nanmean(KAa{n}),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    for e=1:size(KAp{n},1)
        isn = ~isnan(KAa{n}(e,:));
        X = command{n}(e,isn);
        Y = KAa{n}(e,isn);
        KAap{n}(e,:) = lsqcurvefit(bfun,bp0,X(2:end),Y(2:end),[-inf,-inf],[inf,inf],opts);
%         plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%         scatter(command{n}(e,:)+n,KAa{n}(e,:),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    aplt(n) = plot(-80:20,bfun(mean(KAap{n}),-80:20),'Color',colorss(n,:));hold on
end

axs0.Title.String = 'Activation Curve';
axs0.YLabel.String = 'Activation';
axs0.XLabel.String = 'Command potential (mV)';
axs0.XLim = [-80,20];
axs0.YLim = [-0.1 1];
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
%         Yr = prctile(Yp,[5 95]);
%         plot([n , n],Yr,'Color',colorss(n,:));
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
    axes(axs1(3))
%     YaA = [YaA; KAp{n}(:,end,1)./pass{n}(:,6)];
%     scatter(ones(size(KAp{n},1),1)*n, KAp{n}(:,end,1)./pass{n}(:,6),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
%     Yp = prctile(KAp{n}(:,end,1)./pass{n}(:,6),[25 50 75]);
%     rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
%     plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
%     Ys = KAp{n}(:,end,1)/(command{n}(1,end) +70);%./pass{n}(:,6);
    Ys = KAp{n}(:,end,1)/(command{n}(1,end) +70);%./pass{n}(:,6);
    YaA = [YaA; Ys];
    scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Ys,[25 50 75]);
%     Yr = prctile(Ys,[5 95]);
%     plot([n , n],Yr,'Color',colorss(n,:));
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


% scatter time constant
% 
% bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2)));
% bp0 = [-10, 1];
% KAap = cell(1,3);
% aplt = gobjects(1,3);
% for n=1:3
%     Y = mean(KAp{n}(:,:,2));
%     iqr = prctile(KAa{n},[25 75]);
%     err = nanstd(KAa{n})/sqrt(size(KAa{n},1));
%     err = [-err;err] + nanmean(KAa{n});
% 
%     Xx = command{n}(1,:)+n-2;
%     scatter(Xx,nanmean(KAa{n}),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
%     plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
%     for e=1:size(KAp{n},1)
%         isn = ~isnan(KAa{n}(e,:));
%         X = command{n}(e,isn);
%         Y = KAa{n}(e,isn);
%         KAap{n}(e,:) = lsqcurvefit(bfun,bp0,X(2:end),Y(2:end),[-inf,-inf],[inf,inf],opts);on
%     end
%     aplt(n) = plot(-80:20,bfun(mean(KAap{n}),-80:20),'Color',colorss(n,:));hold on
% end

% axs0.Title.String = 'Activation Curve';
% axs0.YLabel.String = 'Activation';
% axs0.XLabel.String = 'Command potential (mV)';
% axs0.XLim = [-80,20];
% axs0.YLim = [-0.1 1];
% axs0.TickDir = 'out';





%time constant averaged across voltages

itc = [0.035, 0.035, 0.14];

axs2 = axes('Position',[0.85 0.55 0.12 0.4]);
colorss = [0.5,0,0;0,0.5,0;0,0,0.5];
Ya = zeros(0,1);
Yg = cell(0,1);
for n=1:3
    Y = median(KAp{n}(:,:,2),2);
    Ya = [Ya;Y]; %#ok<AGROW>
    Yg = [Yg; repmat({names{n}},length(Y),1)];%#ok<AGROW>
    scatter(ones(size(KAp{n},1),1)*n + (rand(size(KAp{n},1),1)-0.5)*0.1,Y,20,'d',...
        'MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Y(:),[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));hold on
    scatter(n,itc(n),'+','MarkerEdgeColor',colorss(n,:));hold on
end
axs2.YLabel.String = 'Time constant (s)';
axs2.XTick = 1:3;
axs2.XTickLabel = ["B51","B64","B8"];
axs2.XLim = [0.2 3.8];
axs2.YLim(1) = 0;

disp('Time constant')
[cta{1},paa(1)] = quickanova(Ya,Yg);
makesig(cta{1},paa(1),axs2);


LaTeX_Expr = 'f(x) = A_1(1-e^{\frac{-x}{\tau_1}})(e^{\frac{-x}{\tau_2}})+b';
text(0, -0.5,['$$' LaTeX_Expr '$$'], 'FontSize',15, 'Color','k', ...
    'HorizontalAlignment','Center', 'VerticalAlignment','Middle', ...
    'Interpreter','Latex');

% %%
% 
% % time constant scatter
% axs21 = axes('Position',[0.05 0.05 0.4 0.4]);
% for n=1:3
%     Y = KAp{n}(:,:,2);
%     scatter(command{n}(:)+n,Y(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
% end
% legend(names)
% axs21.Title.String = 'Inactivation time constant';
% axs21.YLim = [0 0.2];
%

% time constant box
% axs2 = axes('Position',[0.5 0.05 0.2 0.4]);
% for n=1:3
%     color = zeros(14,3);
%     color(:,n) = linspace(0.3,0.9,14);
%     for t=3:size(KAp{n},2)
%         Y = KAp{n}(:,t,2);
%         scatter(ones(size(KAp{n},1),1)*n + (rand(size(KAp{n},1),1)-0.5)*0.1,Y,20,'d',...
%             'MarkerFaceColor',color(t,:),'MarkerEdgeColor',color(t,:));hold on
%     end
%     Y = KAp{n}(:,:,2);
%     Yp = prctile(Y(:),[25 50 75]);
%     rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
%     plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
%     
% end
% axs2.YLabel.String = 'Time constant (s)';
% axs2.YLim = [0,0.3];
% axs2.XTick = 1:3;
% axs2.XTickLabel = names;

% inactivation
opts = optimset('Display','off','Algorithm','levenberg-marquardt');
% axs2(1) = axes('Position',[0.05 0.05 0.3 0.4]);
axes(axs0)

midx = 415:455;
bfun = @(p,x) 1./(1+exp((p(1) - x)/p(2)));
bp0 = [-50, -1];
KAiap = cell(1,3);
KAia = cell(1,3);
for n=1:3
    KAia{n} = squeeze(mean(datai{n}(midx,3,:,3,:)))';
    KAia{n} = KAia{n} - repmat(min(KAia{n},[],2),1,size(KAia{n},2));
    KAia{n} = KAia{n}./repmat(max(KAia{n},[],2),1,size(KAia{n},2));
    %     scatter(command{n}(:)+n-2,KAa{n}(:),20,'d','MarkerFaceColor',colors(n,:),'MarkerEdgeColor',colors(n,:));hold on
    
    iqr = prctile(KAia{n},[25 75]);
    err = nanstd(KAia{n})/sqrt(size(KAia{n},1));
    err = [-err;err] + nanmean(KAia{n});

    Xx = commandi{n}(1,:)+n-2;
    plot(repmat(Xx,2,1),err,'Color',colorss(n,:));hold on
    scatter(Xx,nanmean(KAia{n}),20,'d','MarkerFaceColor','w','MarkerEdgeColor',colorss(n,:));hold on
    for e=1:size(datai{n},5)
        if n==3 && e==5
            KAiap{n}(e,:) = lsqcurvefit(bfun,bp0,commandi{n}(e,1:2:end),KAia{n}(e,1:2:end),[-inf,-inf],[inf,inf],opts);
        else
            KAiap{n}(e,:) = lsqcurvefit(bfun,bp0,commandi{n}(e,:),KAia{n}(e,:),[-inf,-inf],[inf,inf],opts);
        end
%         plot(-60:20,bfun(KAap{n}(e,:),-60:20),'Color',colors(e,:));hold on
%         scatter(command{n}(e,:)+n,KAa{n}(e,:),20,'d','MarkerFaceColor',colors(e,:),'MarkerEdgeColor',colors(e,:));hold on
    end
    plot(-80:20,bfun(mean(KAiap{n}),-80:20),'LineStyle','--','Color',colorss(n,:));hold on
end
% legend(aplt,names,'Location','southeast')
% axs2(1).Title.String = 'Inactivation Curve';
% axs2(1).YLabel.String = 'Activation';
% axs2(1).XLabel.String = 'Command potential (mV)';
% axs2(1).XLim = [-80,20];
legend(aplt,names,'Location','southeast')

%boxplots
axs2(2) = axes('Position',[0.4 0.05  0.12 0.4]);
axs2(3) = axes('Position',[0.55 0.05 0.12 0.4]);

Yiah = zeros(0,2);
Yigh = cell(0,1);
for n=1:3
    Yiah = [Yiah;KAiap{n}]; %#ok<AGROW>
    Yigh = [Yigh;repmat({names{n}},size(KAiap{n},1),1)];%#ok<AGROW>
    for p=1:2
        axes(axs2(p+1))
        scatter(ones(size(KAiap{n},1),1)*n,KAiap{n}(:,p),20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
        Yp = prctile(KAiap{n}(:,p),[25 50 75]);
%         Yr = prctile(Yp,[5 95]);
%         plot([n , n],Yr,'Color',colorss(n,:));
        rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
        plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:))
    end
end 
axs2(2).YLabel.String = 'Half inactivation (mV)';
axs2(3).YLabel.String = 'Slope';
axs2(3).YLim(2) = 0;
set(axs2(2:3),'XTick',1:3);
set(axs2(2:3),'XTickLabel',names);
set(axs2(2:3),'XLim',[0.2 3.8]);
set(axs2,'TickDir','out');
disp('===========Inactivation==============')
disp('Half inactivation')
[cta{1},paa(1)] = quickanova(Yiah(:,1),Yigh);
makesig(cta{1},paa(1),axs2(2))
disp('Slope')
[cta{1},paa(1)] = quickanova(Yiah(:,2),Yigh);
makesig(cta{1},paa(1),axs2(3))



%%
figure('Color','w')
bcolors = makecolor(0.5);
for n=3:-1:1
    fW2 = 4000:10000;
    X = linspace(0,length(fW2)/1e4,length(fW2));
    Y = squeeze(data{n}(fW2,1,end,3,:));
    Y = Y./repmat(max(Y(240:5000,:)),length(fW2),1);
    
    fx = [X fliplr(X)];
    err = std(Y,[],2)/sqrt(size(Y,2));
    err = [-err';err'] + mean(Y');
    fy = [err(1,:) fliplr(err(2,:))];
    fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(n,:));hold on

%     plot(Y,'color',bcolors(p,:));hold on
    plot(X,mean(Y,2),'color',dcolors(n,:));hold on
end
ax = gca;
ax.YLim = [-0.10 1.05];
ax.XLim = [0 0.3];
ax.Box = 'off';

