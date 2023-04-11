%%
% preceeded by KCa_new4MTEA
%%

folder = '\\uthsch-nas.uthouston.edu\ms_nba\ByrneLab\cneveu\Data\2019\';
if strcmp(getenv('COMPUTERNAME'),'MSNBA160674') || strcmp(getenv('COMPUTERNAME'),'MSNBA-C000045')
%     folder = 'C:\Users\cneveu\Desktop\Data\2019\';
%     folder = 'D:\Data\2019\';
%     folder = 'F:\Data\2019\';
%     folder = 'C:\Users\cneveu\Documents\Data\2019\';
    folder = 'C:\Users\cneveu\My Drive\Data\2019\';
end
% folder = 'E:\external_drive\Data\2019\';

neuron = ["B51","B64","B8"];
info=cell(1,3);
data = cell(1,3);
datac = cell(1,3);
logg = '01';
for n=1:length(neuron)
    [~,~,info{n}] = xlsread([folder 'KCa_new_4mTEA.xlsx'],neuron{n});
    info{n} = info{n}(3:end,2:end);
    info{n} = cellfun(@num2str,info{n},'UniformOutput',false);

    [data{n},si] = readsdata(folder,info{n}(:,1:2),1:1.1e4,4.5e3);

    data{n}(:,[1 3],:,3,:) = lowpassf(data{n}(:,[1 3],:,3,:),'Fpass',700,'Fstop',900,'Fs',1e6/si);% 500 700
    
    for e=1:size(info{n},1)
        [datagf,sig] = abfload([folder info{n}{e,3} '.abf'],'verbose',false);
        datagf = datagf(round(str2double(info{n}{e,4})*1e3/sig):round(str2double(info{n}{e,5})*1e3/sig),:);
        pulse = datagf(50:end,4) - datagf(1:end-49,4);
        pup = logg((pulse>0.055)+1);
        pt = strfind(pup,repelem('01',47));
        pdn = logg((pulse<-0.055)+1);
        pe = strfind(pdn,repelem('01',47));
        long = arrayfun(@(x) any(abs(pt-x-25000)<100),pe);       
        pe = pe(long);
        datac{n}(:,:,e) = datagf(pe(1)+50:pe(1)+5000,:);
    end
end
disp('>>>> finished reading data <<<<<')


%% subract NiCd

% runn KCa_new4mMTEA_NiCd
ratio = cell(1,3);
NiCd = cell(1,3);
Control = cell(1,3);
NiCdn = cell(1,3);
Y = cell(1,3);
nW = 4470:5470;
for n=1:3
    NiCd{n} = squeeze(mean(NiCd_all{n}(:,1,:,3,:),5));
    Control{n} = squeeze(mean(data{n}(:,1,:,3,:),5));
    Y{n} = mean(NiCd_all{n}(:,1,:,3,:),5);
    NiCdn{n} = Y{n}./repmat(max(abs(Y{n}(nW,:,:,:,:))),size(NiCd{n},1),1,1,1,1);
    NiCdn{n} = repmat(NiCdn{n},1,1,1,1,size(data{n},5));
    ratio{n} = max(abs(mean(NiCd_all{n}(nW,1,:,3,:),5)))./max(abs(mean(data{n}(nW,1,:,3,:),5)));
    ratio{n}(ratio{n}>1) = 1;
    NiCdn{n} = NiCdn{n}.*repmat(max(abs(data{n}(nW,1,:,3,:))),size(NiCdn{n},1),1,1,1,1).*repmat(ratio{n},size(NiCdn{n},1),1,1,1,1);
end



%% check normalization
figure;n=3;for t=1:9;plot(NiCdn{n}(4000:10000,1,t));hold on;end
%% normalized NiCd

t=9;
for n=1:3
    figure
    for e=1:6
        ax = subplot(2,3,e);
        plot([squeeze(data{n}(4000:10000,1,t,3,e)) squeeze(NiCdn{n}(4000:10000,1,t,1,e)) ]);
        ax.YLim=[-5,300];
        ax.Title.String=neuron{n};
        legend('Control','NiCd');
    end
end


%% current traces

command = cell(1,3);
for n=1:3
    command{n} = round(squeeze(data{n}(5000,2,:,3,:))/10)'*10;
end


sdata = cell(1,3);
for n=1:3
    sdata{n} = data{n}(:,1,:,3,:) - NiCdn{n}(:,1,:,1,:);
end

save('Kca.mat','sdata','command')

for n=1:length(neuron)
    close(findobj(0,'Name', [neuron{n} ' subtractions']))% 
    fig = figure('Name',[neuron{n} ' subtractions'],'NumberTitle','off');
    for e=1:size(info{n},1)
        ax = subplot(2,4,e);
%         plot(400:0.1:1000,squeeze(mean(data{n}(4000:10000,1,:,3,:),5))) 
        Y = squeeze(data{n}(4000:10000,1,:,3,e)) - squeeze(NiCdn{n}(4000:10000,1,:,1,e));
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
dcolors = makecolor(-0.2);
bcolors = makecolor(0.5);

close(findobj(0,'Name', 'subtractions'))% 
fig = figure('Name','subtractions','NumberTitle','off');
for n=1:length(neuron)
    ax = subplot(1,3,n);
    for t=4:size(sdata{n},3)
        Y = squeeze(sdata{n}(4000:10000,1,t,1,:));
        X = 1:size(Y,1);
        mY = mean(Y,2);
        
        fx = [X fliplr(X)];
        err = nanstd(Y,[],2)/sqrt(size(Y,2));
        err = [-err';err'] + nanmean(Y');
        fy = [err(1,:) fliplr(err(2,:))];
        fill(fx,fy,colors(n,:),'EdgeColor','none','FaceColor',bcolors(t-3,:));hold on
        
        plot(mY,'Color',dcolors(t-3,:));hold on 
    end
    
        ax.YLim = [-10, 180];
%     ax.YLim = [-90, 150];
%     ax.Title.String = info{n}{e,1};%['Pre: ', num2str(info{e,3}), ', Post: ', num2str(info{e,4}),', Pre: ', num2str(info{e,5})];
end


%%

figure('Name','scatter','NumberTitle','off','Position',[100 100 300 400]);
names = ["B51","B64","B8"];

%boxplots
axs1 = axes('Position',[0.2 0.1  0.5 0.8]);
Ya = zeros(0,2);
Yg = cell(0,1);
YaA = zeros(0,2);
for n=1:3
    Ys = squeeze(mean(sdata{n}(9200:9301,1,6,1,:)));
%     Ys(Ys<0) = 0;
    Ya = [Ya;Ys]; %#ok<AGROW>
    Yg = [Yg;repmat(names(n),length(Ys),1)];%#ok<AGROW>
    scatter(ones(size(Ys))*n, Ys,20,'d','MarkerFaceColor',colorss(n,:),'MarkerEdgeColor',colorss(n,:));hold on
    Yp = prctile(Ys,[25 50 75]);
    rectangle('Position',[n - 0.25 , Yp(1) , 0.5, Yp(3) - Yp(1)],'EdgeColor',colorss(n,:)); hold on
    plot([n-0.25 , n+0.25],[Yp(2) Yp(2)],'Color',colorss(n,:));
end 
axs1.YLabel.String = 'Current (nA)';
set(axs1,'TickDir','out');
% axs1(3).YLim(2) = 3.5;
set(axs1,'XTick',1:3);
set(axs1,'XTickLabel',names);
set(axs1,'XLim',[0.2 3.8]);

disp('IKCa Peak current')
[cta{1},paa(1)] = quickanova(Ya(:,1),Yg);
makesig(cta{1},paa(1),axs1(1));
