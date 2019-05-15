%bulk import of sialidosis data using fieldtrip toolbox
%last updated 8/23/2018 - Patrick McGurrin
clear all; clc; close all

%what is the patient ID? - only one at a time
subjNum = 9;


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

load(strcat('tremor_',PTid,'.mat'));

titles = {'preStim 1';'preStim 2';'during Stim';'postStim 1';'postStim 2'};
% organizes out of fieldtrip structure
for posti = 1:2
    for n = 1:size(CleanData,1)
        for triali = 1
            trialData{n,posti}(:,:,triali) = CleanData{n,posti}.trial{1,triali};
        end
    end; 
%     clear n triali
end; 

%clear posti CleanData

%%splits into EEG and EMG data

for posti = 1:2
    for condi = 1:size(trialData,1)
        
        eegData{condi,posti} = trialData{condi,posti}(1:63, 1:30000);
        emgData{condi,posti} = trialData{condi,posti}(64:68,1:30000);
        accData{condi,posti} = trialData{condi,posti}(69:70,1:30000); 
%         clear trialDat
        
    end; 
    
end; 
 
%clear posti trialData condi

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

srate = 1000;

min_frex = 2;
max_frex = 20;
num_frex = 60

frex = linspace(min_frex,max_frex,num_frex);

wtime    = -4:1/srate:4;
half_wav = (length(wtime -1))/2;

nData    = length(emgData{1,1});
nKern    = length(wtime);
nConv    = nData + nKern -1;

% Full with at half max

fwhm = linspace(1.5,0.7,length(frex));

gwin = zeros(length(wtime),length(frex));

for li = 1:length(frex)
   gwin(:,li) = exp( (-4*log(2)*wtime.^2) ./ fwhm(li)^2 );
   
end



clear li
% create family of wavelets

for li = 1: length(frex) 
sine_wave(li,:) = exp ( 1i*2*pi*frex(li).*wtime );
cmw             = sine_wave(li,:).* gwin(:,li)';

tempX           = fft(cmw,nConv);
cmwX(li,:)      = tempX/max(abs(tempX));

end

clear tempX cmw sine_wave gwin fwhm  nKern  wtime min_frex max_frex num_frex li n triali condi

% and now, convolution
%EMG
post  = zeros(5,5,length(frex)); % no weigt
postW = zeros(5,5,length(frex)); % with  weigt

for tp = 1:5  % task
    
    for chin = 1:5
        
        chan1 = abs(emgData{tp,1}(chin,:));
        chan2 = abs(emgData{tp,2}(chin,:));
        
        chan1 = fft(chan1,nConv)/nData;
        chan2 = fft(chan2,nConv)/nData;
        
for     li = 1: length(frex) 
        as1 = ifft(chan1.* cmwX(li,:));
        as1 = as1(half_wav-1:end - half_wav);
        post(tp,chin,li)= mean(abs(as1));
        
        as2 = ifft(chan2.* cmwX(li,:));
        as2 = as2(half_wav-1:end - half_wav);
        postW(tp,chin,li)= mean(abs(as2));
    
end


    end
end
clear as1 as2 chan1 chan2 chin ans foo li tp

%ACC
acc  = zeros(5,2,length(frex));  % no weigt
accW = zeros(5,2,length(frex));  % with weigt

for tp = 1:5  % task
    
    for chin = 1:2
        
        chan1 = fft(accData{tp,1}(chin,:),nConv)/nData;
        chan2 = fft(accData{tp,2}(chin,:),nConv)/nData;
        
for     li = 1: length(frex) 
        as1 = ifft(chan1.* cmwX(li,:));
        as1 = as1(half_wav-1:end - half_wav);
        acc(tp,chin,li)= mean(abs(as1));
        
        as2 = ifft(chan2.* cmwX(li,:));
        as2 = as2(half_wav-1:end - half_wav);
        accW(tp,chin,li)= mean(abs(as2));
    
end


    end
end
clear as1 as2 chan1 chan2 chin ans foo li tp nData


%
figure;
 set(gcf,'color','w');
for ti = 1:5
subplot(221)
plot(frex,squeeze(acc(ti,2,:)),'LineWidth',4)
hold on
legend(titles)
title('left NO weight')
subplot(222)
plot(frex,squeeze(acc(ti,1,:)),'LineWidth',4)
hold on
legend(titles)
title('right NO weight')
subplot(223)
plot(frex,squeeze(accW(ti,2,:)),'LineWidth',4)
hold on
legend(titles)
title('left with weight')
subplot(224)
plot(frex,squeeze(accW(ti,1,:)),'LineWidth',4)
hold on
legend(titles)
title('right with weight')

end
suptitle(['ACC, frex domain, subject ' num2str(subjNum) ])




figure;
set(gcf,'color','w');
for ti = 1:5

subplot(321)
plot(frex,squeeze(acc(ti,2,:)),'LineWidth',4)
hold on
legend(titles)
title('ACC left NO weight')

subplot(322)
plot(frex,squeeze(acc(ti,1,:)),'LineWidth',4)
hold on
legend(titles)
title('ACC right NO weight')

subplot(323)
plot(frex,squeeze(post(ti,4,:)),'LineWidth',4)
hold on
legend(titles)
title('FCR right NO weight')

subplot(324)
plot(frex,squeeze(post(ti,5,:)),'LineWidth',4)
hold on
legend(titles)
title('ECR right NO weight')

subplot(325)
plot(frex,squeeze(post(ti,2,:)),'LineWidth',4)
hold on
legend(titles)
title('FCR left NO weight')

subplot(326)
plot(frex,squeeze(post(ti,3,:)),'LineWidth',4)
hold on
legend(titles)
title('ECR left NO weight')


end
suptitle(['Tremor, frex domain, subject ' num2str(subjNum) ])
saveas(gcf,['sub_' num2str(subjNum) '_tremor.fig' ])
saveas(gcf,['sub_' num2str(subjNum) '_tremor.png' ])

%% Set range in wich amplitud is going to be look for
r1 = 5.5; % in Hz
r2 = 7; % in Hz

[~,r1] = min(abs(frex - r1));
[~,r2] = min(abs(frex - r2));


reACC   = reshape(acc,  10,60);
reACCW  = reshape(accW, 10,60);
rePost  = reshape(post, 25,60);
rePostW = reshape(postW,25,60);

ampACC  = zeros(10,1);
ampACCW = zeros(10,1);
for ti = 1:10
   ampACC(ti)  = max(reACC(ti,r1:r2));
   ampACCW(ti) = max(reACCW(ti,r1:r2));
    
end

ampPost  = zeros(25,1);
ampPostW = zeros(25,1);
for ti = 1:10
   ampPost(ti)  = max(rePost(ti,r1:r2));
   ampPostW(ti) = max(rePostW(ti,r1:r2));
    
end

clear r1 r2 ti reACC reACCW rePost rePostW posti



ampACC   = reshape(ampACC,  5,2);
ampACCW  = reshape(ampACCW, 5,2);
ampPost  = reshape(ampPost, 5,5);
ampPostW = reshape(ampPostW,5,5);


figure;
set(gcf,'color','w');
subplot(121)
bar(ampACC,'DisplayName','ampACC')
set(gca,'xticklabel',titles);
legend({'right'; 'left'})
title('no Weigth')
subplot(122)
bar(ampACCW,'DisplayName','ampACCW')
set(gca,'xticklabel',titles);
legend({'right'; 'left'})
title('With Weigth')
suptitle(['Amplitude tremor(ACC), subject ' num2str(subjNum) ])


%%
posti = 1;
%1 = posture
%2 = posture w/ 1 pound


figure(); set(gcf,'color','w');


frexRange = [2.5 20];


plotNums = [1 2 3;
            4 5 6;
            7 8 9;
            10 11 12;
            13 14 15];

for condi = 1:length(emgData)

    for n = 1:3 %ACC, ECR, FCR
        subplot(5,3,plotNums(condi,n))
        
        %selection of left or right data
       
        if n == 1
            dat = accData{condi,posti}(n+1,1:30000)'; %+1 because first channel is actually APB
        else
            dat = emgData{condi,posti}(n+1,1:30000)'; %+1 because first channel is actually APB
        end
        
        if n == 2 || n == 3
            dat = abs(dat);
        end
        
        hz = linspace(0,sRate/2,floor(length(dat)/2)+1);
        
        
        %processing steps
        %f  = fft(dat(:,1)/length(dat));
        f  = fft(log10(dat(:,1))/length(dat)).^2;
        
        for nn = 1:length(frexRange)
            [x,y] = min(abs(hz-frexRange(1,nn)));
            fr(nn,1) = y;
        end; clear nn y
        
        %now plot
        x = hz';
        %y = abs(f(1:length(hz))*2);
        %y = smooth(abs(f(1:length(hz))*2),'lowess');
        y = abs(f(1:length(hz))*2);
        y = envelope(y,20,'rms');
        
        
        [yy xx] = max(y(fr(1):fr(2),1));
        
        plot(x,y,'b'); hold on;
        %plot(x,smooth(y,'lowess'),'b'); hold on;
        
        if  n == 1
            tremPower(:,condi) = y;
        end
        
        %if nums(n,1) == 1 || nums(n,1) == 2
        scatter(x(ceil(xx)+fr(1),1),yy,'r','MarkerEdgeColor','r','MarkerFaceColor','r');
        text(x(ceil(xx)+fr(1),1),yy,['(', num2str(x(ceil(xx)+fr(1),1)), ', ', num2str(yy), ')'])
        %end
        
        xlim([2 20]);
        %set(gca,'ylim',[-15 0])
        

        if condi == 1 && n == 1
            title('right ACC');
        elseif condi == 1 && n == 2
            title('right FCR');
        elseif condi == 1 && n == 3
            title('right ECR');
        end
        
        if n == 1
            ylabel(titles{condi,1});
        end
        
        set(gca,'box','off')
        
    end; clear n
    
end; clear condi

fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',14,'FontName','Helvetica'); hold on; clear fig dat

%%
t1 = 11000;
t2 = 12000;

figure(); set(gcf,'color','w');

for condi = 1:length(emgData)

    for n = 1:3 %ACC, ECR, FCR
        subplot(5,3,plotNums(condi,n))
        
        %selection of left or right data
        if n == 1
            dat = accData{condi,posti}(n+1,t1:t2)'; %+1 because first channel is actually APB
        else
            dat = emgData{condi,posti}(n+1,t1:t2)'; %+1 because first channel is actually APB
        end
        
        plot(dat,'k'); hold on;
   
       
        if condi == 1 && n == 1
            title('right ACC');
        elseif condi == 1 && n == 2
            title('right FCR');
        elseif condi == 1 && n == 3
            title('right ECR');
        end
        
        if n == 1
            ylabel(titles{condi,1});
        end
        
        set(gca,'box','off')
        
    end; clear n
    
end; clear condi

fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',14,'FontName','Helvetica'); hold on; clear fig dat


%% overlay tremor powers - leaving out the during stim one
figure(); set(gcf,'color','w');
plotNums = [1 2 4 5]';

for condi = 1:length(plotNums)
    plot(x,tremPower(:,plotNums(condi,1))); hold on;
    xlim([2 20]);
end; clear condi

set(gca,'box','off')
fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',14,'FontName','Helvetica'); hold on; clear fig dat

f=get(gca,'Children');
%legend([f(5),f(4),f(3),f(2),f(1)],'bs1','bs2','ds','ps1','ps2','location','northoutside','Orientation','horizontal');clear fig f 
legend([f(4),f(3),f(2),f(1)],'bs1','bs2','ps1','ps2','location','northoutside','Orientation','horizontal');clear fig f 
legend boxoff
    