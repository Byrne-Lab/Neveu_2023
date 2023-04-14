function datac = readpass(folder,info,fig) 
%this function returns the window of data for the capacitance measurements
% fig will generate a figure

% [~,~,info] = xlsread('D:\Data\2019\Rtype.xlsx','B51');
% info = info(2:end,4:6);
% info = cellfun(@num2str,info,'UniformOutput',false);
% datac = readpass('D:\Data\2019\',info,true);

if nargin<3
    fig = false;
end
logg = '01';
% datac = nan(25000,2,size(info,1),25);
datac = nan(1.5e5,2,size(info,1),25);
if fig;figure;end
for e=1:size(info,1)
    [dataf,sig] = abfload([folder info{e,1} '.abf'],'verbose',false);
    
    if strcmp(info{e,1}(1:2),'21')
        idx = [3 4];
    else
        idx = [2 3];
    end
    
    datagf = dataf(round(str2double(info{e,2})*1e3/sig):round(str2double(info{e,3})*1e3/sig),idx);
    pulse = datagf(50:end,2) - datagf(1:end-49,2);
    if strcmp(info{e,1},'21601013')
        pt = [959312;1033578;1052230;1071220;1089871;1107844;1125139;1141756;1158373;1174650;1191267;1206527;1223822] - round(str2double(info{e,2})*1e3/sig);
    elseif strcmp(info{e,1},'21601004')
        pt = [120447;145876;169244;199141;242784;296048;336598;382646;443471;498454;536598;589519];
    else
        pup = logg((pulse>0.055 | pulse<-0.055)+1);
        pt = strfind(pup,repelem('01',47));
        long = arrayfun(@(x) any(abs(pt-x-26000)<100 | abs(pt-x-25000)<100 | abs(pt-x-5000)<100),pt);       
        pt = pt(long);
    end
    
    
    for t=1:length(pt)
%         W = pt(t)+50: pt(t)+25049;
        W = pt(t)+50: pt(t)+1.5e5+49;
        W(W>length(datagf)) = [];
        datac(1:length(W),:,e,t) = datagf(W,:);
    end

    if fig
       subplot(size(info,1),1,e)
       plot(datagf(:,2));hold on
%        plot(diff(datagf(:,2)));hold on
       plot(pulse);hold on
       scatter(pt,ones(size(pt))*mean(datagf(:,2)));hold on
    end
end
datac(:,:,:,all(isnan(datac),1:3)) = [];
end

