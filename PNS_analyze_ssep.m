%bulk import of sialidosis data using fieldtrip toolbox
%last updated 8/23/2018 - Patrick McGurrin
clear all; clc; close all

%what is the patient ID? - only one at a time
subjNum = 6;


if subjNum == 4
    PTid = '4_TS00013';
elseif subjNum == 5
    PTid = '5_TS00025';
elseif subjNum == 6
    PTid = '6_TS00037';
elseif subjNum == 7
    PTid = '7_TS00002';
elseif subjNum == 9
    PTid = '9_TS00010';
end

%everything else below will run on its own!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chan2use = {'fc3';'c3';'cp3';'p3'}; %for multi EEG plot
chan2useSingle = {'cp3'}; %for SSEP c wave plot - single EMG/EEG

%%load data

%where is the data? -- and call in relevant folders for data/toolbox access
if ismac == 1
    data_Loc = strcat('/Volumes/shares/DIRFS1/Protocol 17-N-0035/01_PNS Substudy/Data/',PTid,'/');
else
    data_Loc = strcat('\\nindsdirfs\Shares\DIRFS1\Protocol 17-N-0035\01_PNS Substudy\Data\',PTid,'\');
end
cd(data_Loc)

load(strcat('SSEP_',PTid,'.mat'));

titles = {'preStim 1';'preStim 2';'postStim 1';'postStim 2'};


% organizes out of fieldtrip structure
for n = 1:size(CleanData,1)
    for triali = 1:nTrials
        trialData{n,1}(:,:,triali) = CleanData{n,1}.trial{1,triali};
    end
end; clear n triali CleanData 

%%splits into EEG and EMG data
for condi = 1:size(trialData,1)

    eegData{condi,1} = trialData{condi,1}(1:63,:,:);
    emgData{condi,1} = trialData{condi,1}(64,:,:); clear trialDat

end; clear trialData condi

%modify some values for plotting purposes below
winSizeBefore = winSizeBefore*1000;
winSizeAfter = winSizeAfter*1000;

chan =find(ismember(eegNames,'CP3'));

figure;
for ti = 1:4
subplot(4,1,ti)
dat = squeeze(eegData{ti,1}(chan,:,:));
mesh(dat);

end

%% Correct and reject trials

% for subjet 7
eegData{3,1} = eegData{3,1}(:,:,200:end);

%%
% Time for the cursors
val1 = 18;
val2 = 27;

%determine data for plotting
for n = 1:length(eegNames)
    chan(n,1) = strcmpi(chan2useSingle,eegNames{n,1});
end; clear n; p = find(chan(:,1) == 1); clear chan

figure(); set(gcf,'color','w'); hold on;


ttime = -50:199;

plotNum = [1;3;5;7];

for ploti = 1:4
    subplot(4,2,plotNum(ploti))
    
    chan1 = squeeze(eegData{ploti,1}(p,:,:));
    
    for numi = 1:size(chan1,2)
        chan1(:,numi) = detrend(chan1(:,numi));
        chan1(:,numi) = chan1(:,numi) - mean(chan1(1:25,numi));
        
        
        
        
    end; clear numi
    
    chan1  = mean(chan1,2);
    newEEG(ploti,:) = mean(chan1,2);
    
    ssepAmp(ploti) = peak2peak(chan1(val1+50:val2+50,1));
    
    
    plot(ttime,chan1,'k'); hold on;
    title(titles{ploti,1})
    %analysis lines
    line([val1 val1],[min(chan1) max(chan1)],'color',[.5 .5 .5]);hold on;
    line([val2 val2],[min(chan1) max(chan1)],'color',[.5 .5 .5]);hold on;

    line([-50 200], [0 0],'color',[.5 .5 .5]);hold on;
    set(gca,'Ydir','reverse')
    
    xlim([-50 200]);
     
    ylabel(chan2useSingle);
    set(gca,'box','off');
    
    if ploti == 4
        xlabel('Time (ms)');
    else
        set(gca,'xcolor','w');
    end
     ylim([max(ssepAmp)*-2 max(ssepAmp)*2])
%     xticks([0:25:winSizeAfter+winSizeBefore])
%     xticklabels({'-50','-25','0','25','50','75','100','125','150','175','200','225','250'})
    suptitle (['sub ' num2str(subjNum) ' SSEP' ])
    
    
end; clear ploti

fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',12,'FontName','Helvetica'); hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])

ssepAmp

subplot(4,2,[2;4;6;8])
for pli = 1:4
    
   plot(ttime,newEEG(pli,:).*-1)
   hold on
end
legend(titles)

saveas(gcf,['sub_' num2str(subjNum) '_SSEP.fig' ])
saveas(gcf,['sub_' num2str(subjNum) '_SSEP.png' ])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val1 = 68;
val2 = 80;

%determine data for plotting
for n = 1:length(eegNames)
    chan(n,1) = strcmpi(chan2useSingle,eegNames{n,1});
end; clear n; p = find(chan(:,1) == 1); clear chan

figure(); set(gcf,'color','w'); hold on;

%plot APB
plotNum = [1 3 5 7]';
for ploti = 1:4
    subplot(4,2,plotNum(ploti,1))
    
    x = squeeze(emgData{ploti,1}(1,:,:));
    
    for numi = 1:size(x,2)
        x(:,numi) = detrend(x(:,numi));
        x(:,numi) = x(:,numi) - mean(x(1:25,numi));
    end; clear numi
    
    x = (mean(x,2));
    
    plot(x,'b'); hold on;
    
    line([0 winSizeAfter+winSizeBefore], [0 0],'color',[.5 .5 .5]);hold on;
    
    xlim([0 winSizeAfter+winSizeBefore]);
    set(gca,'box','off'); %,'xcolor','w');
    ylabel(titles{ploti,1})
    ylim([-50 50]);
    
    if ploti == 4
        xlabel('Time (ms)');
    else
        set(gca,'xcolor','w');
    end
%     
%     xticks([0:25:winSizeAfter+winSizeBefore])
%     xticklabels({'-50','-25','0','25','50','75','100','125','150','175','200','225','250'})
%     
end

%plot eeg channel
plotNum = [2 4 6 8]';
for ploti = 1:4
    subplot(4,2,plotNum(ploti,1))
    
    chan1 = squeeze(eegData{ploti,1}(p,:,:));
    
    for numi = 1:size(chan1,2)
        chan1(:,numi) = detrend(chan1(:,numi));
        chan1(:,numi) = chan1(:,numi) - mean(chan1(1:25,numi));
        
        %ssepAmp(numi,ploti) = rms(chan1(val1:val2,numi));
        ssepAmp(numi,ploti) = trapz(abs(chan1(val1:val2,numi)));
        
    end; clear numi
    
    chan1 = mean(chan1,2);
    
    plot(chan1,'k'); hold on;
    
    %analysis lines
    line([val1 val1],[min(chan1) max(chan1)],'color',[.5 .5 .5]);hold on;
    line([val2 val2],[min(chan1) max(chan1)],'color',[.5 .5 .5]);hold on;

    line([0 winSizeAfter+winSizeBefore], [0 0],'color',[.5 .5 .5]);hold on;
    set(gca,'Ydir','reverse')
    
    xlim([0 winSizeAfter+winSizeBefore]);
    ylabel(chan2useSingle);
    set(gca,'box','off');
    
    if ploti == 4
        xlabel('Time (ms)');
    else
        set(gca,'xcolor','w');
    end
%     
%     xticks([0:25:winSizeAfter+winSizeBefore])
%     xticklabels({'-50','-25','0','25','50','75','100','125','150','175','200','225','250'})
    
end; clear ploti

fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',12,'FontName','Helvetica'); hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])

meanAmp = mean(ssepAmp)