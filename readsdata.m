function [adata,si,commando] = readsdata(folder,info,W,mv,div,chan,filter)
% W is the window of data to return in vector format example units of datapoints : 1:1e4
% info is the cell array from the excel spreadsheet
% mv is the index timepoint to measure voltage 
% div is the division increment of segregating voltage steps 10 == 10mV
% increments
% chan = the channels to be recorded [Im mV]  optional input
% filter = do you want the data be filtered prior to subtracting true/false

if nargin<5
    div=10;
end

if nargin<6
    chan = 1:3;
end

if nargin<7
    filter=false;
end

tic
command = -140:div:60;
adata = nan(length(W),length(chan),length(command),3,size(info,1));
for e=1:size(info,1)
    disp([num2str(e) '  ' info{e,1} '  ' 'Time: ' num2str(toc) ' s'])
    [data1,si] = abfload([folder info{e,1} '.abf'],'verbose',false);
    mvc = round(squeeze(data1(mv,chan(2),:))/div)*div;
    mvidx = arrayfun(@(x) find(command==x),mvc);
    adata(:,:,mvidx,1,e) = data1(W,chan,:);
    if size(info,2)==1
        adata(:,:,mvidx,2,e) = 0;
        adata(:,2,mvidx,2,e) = adata(:,chan(2),mvidx,1,e);
    elseif ~strcmpi(info{e,2},'nan')
        [data2,~ ] = abfload([folder info{e,2} '.abf'],'verbose',false);
        adata(:,:,mvidx,2,e) = data2(W,chan,:);
    end   
end

if filter
    adata(:,[1 3],:,1:2,:) = lowpassf(adata(:,[1 3],:,1:2,:),'Fpass',100,'Fstop',250,'Fs',1e6/si);  
end

adata(:,:,:,3,:) = adata(:,:,:,1,:) - adata(:,:,:,2,:);
adata(:,2,:,3,:) = adata(:,2,:,1,:);
adata(:,:,all(isnan(adata(1,1,:,1,:)),5),:,:) = [];
commando = (round(squeeze(adata(mv-W(1)+1,2,:,1,:))/div)*div)';
end