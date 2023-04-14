function checksim(ch)

folder = 'C:\Users\cneveu\Documents\GitHub\CPG2020\smu';% change this to the location of the simulation files

% uncomment the channel you want to compare to the model. Note, you may
% have to change the file name of the emprical data in the plotmodel
% function below.
if nargin==0
    ch = strings(1,0);
%     ch = [ch, "K"];
%     ch = [ch, "KA"];
%     ch = [ch, "CaR"];
%     ch = [ch, "CaL"];
%     ch = [ch, "Kca"];
%     ch = [ch, "KTEA"];
%     ch = [ch, "Napp"];
    ch = [ch, "HCN"];
end

if any(ch=="K")
    chan = 'K';
    isi = 60;
    fW2 = 4000:16000;
    xlim = [0 1000];
    ylim = [-10 300];
    xylim2 = [0 1000;0 120];
    tpos = 50;
    offset = 5;
    baseline = true;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,xylim2,folder)
    % set(gcf,'Position',[0 200 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="KA") 
    chan = 'KA';
    isi = 4;
    xlim = [0 150];
    ylim = [-10 Inf];
    tpos = 50;
    offset = -5;
    fW2 = 4000:10000;
    baseline = true;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
    % set(gcf,'Position',[0 200 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="CaR") 
    chan = 'CaR';
    isi = 30;
    offset = 1;
    fW2 = 4000:10000;
    xlim = [0 280];%275
    tpos = 50;
    ylim = [-15 5];
    baseline = false;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
    % set(gcf,'Position',[50 150 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="CaL") 
    chan = 'CaL';
    isi = 30;
    offset = 1;
    fW2 = 4000:10000;
    xlim = [0 160];
    tpos = 50;
    ylim = [-10 5];
    baseline = true;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
%     set(gcf,'Position',[100 100 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="Kca") 
    chan = 'Kca';
    isi = 30;
    tpos = 50;
    offset = 2;
    fW2 = 4000:10000;
    xlim = [0 600];
    ylim = [-5  140];
    baseline = true;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
    set(gcf,'Position',[100 120 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="KTEA") 
    chan = 'KTEA';
    isi = 30;
    tpos = 50;
    offset = 2;
    fW2 = 4000:10000;
    xlim = [0 600];
    ylim = [-5  125];
    baseline = true;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
    set(gcf,'Position',[100 120 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="Napp") 
    chan = 'Napp';
    isi = 30;
    offset = 7;
    tpos = 50;
    fW2 = 3500:13500;
    xlim = [0 1000];
    ylim = [-15  5];
    baseline = false;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
%     set(gcf,'Position',[100 120 1600 800])
    set(gcf,'Position',[200 200 300 500])
end

if any(ch=="HCN") 
    chan = 'HCN';
    isi = 80;
    fW2 = 25000:65000;
    offset = 752;
    tpos = 3000;
    xlim = [0 4000];
    ylim = [-7  2];
    baseline = false;
    plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,[],folder)
    set(gcf,'Position',[100 120 1600 800])
end
end


function plotdata(chan,isi,xlim,ylim,offset,baseline,tpos,fW2,xylim2,folder)

close(findobj(0,'Name', chan))% 
figure('Name',chan,'NumberTitle','off','Color','w');


if strcmp(chan,'Kca')
    load([chan '.mat'],'sdata','command')
    data = sdata;
else
    load([chan '.mat'],'data','command')
end
%     for n=1:3
%         pdata = data{n};
%         data{n} = nan([size(pdata,1,2,3) 3 size(pdata,4)]);
%         data{n}(:,:,:,3,:) = pdata;
%     end



if strcmp(chan,'CaR')
    fname = [folder '\' chan '1.smu.out'];
    sdata = dlmread(fname);
    fname = [folder '\' chan '2.smu.out'];
    sdata2 = dlmread(fname);
    sdata(:,[2 4 6]) = sdata(:,[2 4 6]) - sdata2(:,[2 4 6]);
elseif strcmp(chan,'Kca')
    fname = [folder '\' chan '.smu.out'];
    sdata = dlmread(fname);
%     sdata(:,[2 4 6]) = sdata(:,[2 4 6]) + sdata(:,8:10);
else
    fname = [folder '\' chan '.smu.out'];
    if exist(fname,'file')
        sdata = dlmread(fname);
    else
        warning(['The output file for ' fname ' not found'])
        return
    end
end
str = fileread([fname '.head']);
hout = string(regexp(str,'(Ivd|V|IVDxREG)[\w+.\w*','match')');
chidx = contains(hout,chan);


nname = string(regexp(str,'(?<=V[)\w+','match')')';
xx = sdata(:,1);
si = diff(xx(1:2));


if strcmp(chan,'Kca')
    nn = [1:3];
elseif strcmp(chan,'HCN')
    nn = 3;
else
    nn = 1:3;
end
X = (0:length(fW2)-1)/10;
colors = makecolor;
for n=nn
    nam = regexprep(nname{n},'[sa]','');
    nidx = [false; contains(hout,nam) & chidx];
    nvidx = [false; contains(hout,nam) & ~chidx];
    sdata2 = sum(sdata(:,nidx),2);
    sfW2 = round(fW2(1)*1e-4/si:fW2(end)*1e-4/si)+offset;
    xs = xx(1:length(sfW2))*1e3;
    subplot(1,3,n)
    for p=1:size(data{n},3)
        if matches(chan,("K"|"HCN"|"Kca")) 
            Y = squeeze(mean(data{n}(fW2,1,p,1,:),5,'omitnan'));
        else
            Y = squeeze(mean(data{n}(fW2,1,p,3,:),5,'omitnan'));
        end
        
        if baseline
            Y = Y - mean(Y(1:10,:),'omitnan');
        end
        Ys = sdata2(sfW2 + (p-1)*isi/si - round(100*1e-4/si));
        mv = sdata(sfW2(1)+ (p-1)*isi/si + round(1000*1e-4/si),nvidx);
        if ~(n==1 && mv>0 && strcmp(chan,'K_A'))
            plot(X,Y,'Color',colors(p,:));hold on
            plot(xs, Ys,'Color',colors(p,:),'LineStyle','--');hold on% Kx+KCa
            xp = round(length(Ys)*0.75);
%             if strcmp(chan,'K')
%                 text(831,Ys(find(xs>831,1))+5,num2str(mv),'Color',colors(p,:),'VerticalAlignment','middle');hold on
%             elseif strcmp(chan,'Kca')
%                 text(400,Ys(find(xs>400,1))+1,num2str(mv),'Color',colors(p,:),'VerticalAlignment','middle');hold on
%             else
            text(tpos,Ys(find(xs>tpos,1))+0.3,num2str(mv),'Color',colors(p,:),'VerticalAlignment','middle');hold on
%             end
        end
    end
    ax = gca;
    ax.Title.String = nname{n};
    ax.YLabel.String = 'Current (nA)';
    ax.XLabel.String = 'Time (ms)';
    ax.Box = 'off';
    ax.TickDir = 'out';
    if ~isempty(ylim)
        ax.YLim = ylim;
    end
    if ~isempty(xlim)
        ax.XLim = xlim;
    end
end

close(findobj(0,'Name', [chan 'compare']))% 
figure('Name',[chan 'compare'],'NumberTitle','off','Color','w');

p = 4;
for n=nn
    if matches(chan,("K"|"HCN"|"Kca")) 
        Y = squeeze(data{n}(fW2,1,p,1,:));
    else
        Y = squeeze(data{n}(fW2,1,p,3,:));
    end
    
    fx = [X fliplr(X)];
    err = std(Y,0,2,'omitnan')/sqrt(size(Y,2));
    err = [-err';err'] + mean(Y,2,'omitnan')';
    fy = [err(1,:) fliplr(err(2,:))];

    fill(fx,fy,'r','EdgeColor','none','FaceColor',colors(n,:),'FaceAlpha',0.5);hold on          
    plot(X,nanmean(Y,2)','Color',colors(n,:),'LineWidth',3);hold on 
end

nam = regexprep(nname{n},'[sa]','');
nvidx = [false; contains(hout,nam) & ~chidx];
mv = sdata(sfW2(1)+ (p-1)*isi/si + round(1000*1e-4/si),nvidx);
ax = gca;
ax.Title.String = ['Current response at ' num2str(mv)];
ax.YLabel.String = 'Membrane Current (nA)';
ax.XLabel.String = 'Time (ms)';
if nargin>8 && ~isempty(xylim2)
    ax.YLim = xylim2(2,:);
    ax.XLim = xylim2(1,:);
end

end
