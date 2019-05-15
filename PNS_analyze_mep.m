%bulk import of sialidosis data using fieldtrip toolbox
%last updated 8/23/2018 - Patrick McGurrin
clear all; clc; close all

%what is the patient ID? - only one at a time
subjNum = 9; %no meps for subj 4 (from bv - but on spike) 


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
%%load data

%where is the data? -- and call in relevant folders for data/toolbox access
if ismac == 1
    data_Loc = strcat('/Volumes/shares/DIRFS1/Protocol 17-N-0035/01_PNS Substudy/Data/',PTid,'/');
else
    data_Loc = strcat('\\nindsdirfs\Shares\DIRFS1\Protocol 17-N-0035\01_PNS Substudy\Data\',PTid,'\');
end
cd(data_Loc)

load(strcat('mep_',PTid,'.mat'));

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val1 = 20;
val2 = 38;
ttime = -50:199;

figure(); set(gcf,'color','w'); hold on;

%plot APB
for ploti = 1:4
    subplot(4,1,ploti)
    
    x = squeeze(emgData{ploti,1}(1,:,:));
    
    for numi = 1:size(x,2)
        x(:,numi) = detrend(x(:,numi));
        x(:,numi) = x(:,numi) - mean(x(1:25,numi));
        
        %mepAmp(numi,ploti) = rms(x(val1:val2,numi));
        
        
    end; clear numi
    
    x = (mean(x,2));
    Amp(ploti) = peak2peak(x(val1+50:val2+50,1));
    
    plot(ttime,x,'b'); hold on;
    
    %analysis lines
    line([val1 val1],[min(x) max(x)],'color',[.5 .5 .5]);hold on;
    line([val2 val2],[min(x) max(x)],'color',[.5 .5 .5]);hold on;

    
    line([-50 200], [0 0],'color',[.5 .5 .5]);hold on;
    
    xlim([-50 200]);
%     ylim ([-800 800])
    set(gca,'box','off'); %,'xcolor','w');
    ylabel(titles{ploti,1})
    %ylim([-50 50]);
    
    if ploti == 4
        xlabel('Time (ms)');
    else
        set(gca,'xcolor','w');
    end
    
    ylim([max(Amp)*-1 max(Amp)])
%     
%     xticks([0:25:winSizeAfter+winSizeBefore])
%     xticklabels({'-50','-25','0','25','50','75','100','125','150','175','200','225','250'})
%     

suptitle (['sub ' num2str(subjNum) ' MEP' ])
end
fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',12,'FontName','Helvetica'); hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])
saveas(gcf,['sub_' num2str(subjNum) '_MEP.fig' ])
saveas(gcf,['sub_' num2str(subjNum) '_MEP.png' ])


Amp = Amp'