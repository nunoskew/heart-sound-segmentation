function finalmaxs=heartSoundSegmentation2(str,sound,level,fs)%,filename,strfinal[consecutives,avgSim,avgEnergy,avgDifmaxis,avgMinis,avgPikis,avgCvars,avgdifmaxs,avgMaxiwavs,stdMaxiwavs,medianMaxs,stdSim]
% sound: One line of the mtx.mat which is one recorded heartbeat zeropadded
% until it has length which is a power of 2.
% str: which motherwavelet will it be used for the SWT
% at what level will be focusing on to segment the signal by the inflection
% points

close all
z=zeros(1,2^16);
z(1:length(sound))=sound;
z=[z,length(sound)];
% mtx has always lengthof a power of 2 +1 because in the last column of the 
% matrix mtx has the length of the actual signal,so this step ensures 
% that it has a power of 2
a=z(1:end-1);
%z=zeros(1,2^15);
%z(1:length(a))=a;
%a=z;
set(findall(gcf,'type','text'),'FontName','times new roman','fontSize',15)
set(0,'DefaultAxesFontSize',12)
% calls SWT providing the signal, the order of the motherwavelet, and the
% motherwavelet respectively
[swa,swd]=swt(a,10,str);

%store detail coefficients at the desired level
info=swd(level,1:z(end));

% zers stores the changes of the sign of the 
% second derivative of the detail coefficients at the
% desired levels
zers=diff(sign(diff(diff(info))));
% zs stores all the indexes of the non zero values of zers, i.e. the 
% inflection points
zs=find(zers~=0);
% makes sure that the indexes don't surpass the length of the original
% signal, given that it was padded with 0's for the SWT
zs=zs(zs<z(end));
%initializes s that will store the centers of each window delimited by 
%inflection points
s=[];
for(j=1:length(zs)-1)
    s=[s,(zs(j)+zs(j+1))/2];
end

%f=1.7;
%for delimitation purposes we add the 1 and the last point of the original
%signal
if(sum(find(zs==1))==0)
    zs=[1 ,zs];
end
if(sum(find(zs==z(end)))==0)
    zs=[zs,z(end)];
end

bits=cell(1,length(zs)-1);
labels=cell(1,length(zs)-1);
filenames=cell(1,length(zs)-1);
% mtxSim is a similarity matrix of all the segments with each other 
% mtxSim=zeros(length(zs)-1,length(zs)-1);
% %Here we fill the matrix with its similarity values computed either with
% %the energy of the cross correlation between 2 segments or with the dynamic
% %time warping distance between 2 segments (dtw is very slow) 
% for(i=1:length(zs)-1)
%     %i/(length(zs)-1)*100
%     for(j=1:length(zs)-1)
%         if(i<j)
%             mtxSim(i,j)=sum(xcorr(sound(zs(i):zs(i+1)-1),sound(zs(j):zs(j+1)-1)).^2);
%             %mtxSim(i,j)=dtw(sound(zs(i):zs(i+1)-1),sound(zs(j):zs(j+1)-1));
%             mtxSim(j,i)=mtxSim(i,j);
%         end
%     end
% end
% %Normalize by the max value of each line [0,1]
% for(i=1:length(zs)-1)
%     mtxSim(i,:)=mtxSim(i,:)/max(mtxSim(i,:));
% end
% % transform the similarity matrix into a distance matrix 
% mtxSim=1-mtxSim;

% Compute average distance and its standard deviation from the distance 
% matrix
% avgSim=sum(sum(triu(mtxSim)))/(((length(zs)-1)^2)-(length(zs)-1));
% stdSim=sum(std(triu(mtxSim)));

% Initialize metrics: index of the maximum value of each segment,standard
% deviation, maximum value, minimumm value, energy,  coeficient of variation,
% index of the minimum value
maxs=[];
stds=[];
maxis=[];
minis=[];
sumis=[];  
cvars=[];
sumis2=[];
maxiwavs=[];
maxwavs=[];
mins=[];
ant=0;
% fills the metric vectors with their respective values
for(i=1:length(zs)-1)
    plc=find(sound(zs(i):zs(i+1)-1)==max(sound(zs(i):zs(i+1)-1)));
    plc2=find(info(zs(i):zs(i+1)-1)==max(info(zs(i):zs(i+1)-1)));
    plc3=find(info(zs(i):zs(i+1)-1)==min(info(zs(i):zs(i+1)-1)));
    tmpMax=max(sound(zs(i):zs(i+1)-1));
    tmpMaxWav=max(abs(info(zs(i):zs(i+1)-1)));
    tmpMin=min(sound(zs(i):zs(i+1)-1));
    %sumis=[sumis,sum(abs(sound(zs(i):zs(i+1)-1)))/length(sound(zs(i):zs(i+1)-1))];
    cvars=[cvars,std(sound(zs(i):zs(i+1)-1))/mean(sound(zs(i):zs(i+1)-1))];
    sumis2=[sumis2,sum(abs(sound(zs(i):zs(i+1)-1)).^2)];
    maxiwavs=[maxiwavs,tmpMaxWav(1)];    
    maxis=[maxis,tmpMax(1)];
    minis=[minis,tmpMin(1)];
    mins=[mins,plc3+ant];
    maxwavs=[maxwavs,plc2(1)+ant];
    maxs=[maxs,plc(1)+ant];
    ant=zs(i+1);
    %stds=[stds,weightedStd(sound(zs(i):zs(i+1)-1))];
end
%sum of maximum value with the absolute value of the minimum value,
pikis=abs(maxis)+abs(minis);
se=shannonEnergyEnvelope(sound,fs);
hr = HeartRate(se, fs);
nHB=length(info)/(hr*fs)
energies=[];
for(i=1:length(zs)-1)
     energies=[energies,sum(abs(sound(zs(i):zs(i+1)-1)).^2)];
end
energies=energies/max(energies);
nonpeaks=1:length(energies);
Y=pdist([energies(nonpeaks)']);
Z = linkage(Y,'complete');
k2=cluster(Z,2);

pkdtc1=k2==1;
pkdtc2=k2==2;

% finds which cluster has the biggest mean of energies and assumes that is
% the signal cluster's, i.e. which contains the S1's and S2's and the
% clusters, and the other being the noise cluster
if mean(energies(pkdtc1))==max(mean(energies(pkdtc1)),mean(energies(pkdtc2)))
    peaks=pkdtc1;nonpeaks=pkdtc2;
else
    peaks=pkdtc2;nonpeaks=pkdtc1;
end
% sprintf('peaks:%d',length(find(peaks)))
cont=0;
while cont<1%length(find(peaks))<round(nHB)*1.5
    %length(energies(nonpeaks))
    

    ma=energies(nonpeaks);
%     sprintf('length(energies(nonpeaks)):%d',length(energies(nonpeaks)))
    Y=pdist([energies(nonpeaks)']);
    Z = linkage(Y,'average');
    k2=cluster(Z,2);
    pkdtc1=k2==1;
    pkdtc2=k2==2;
    pikis2=ma;
    if mean(pikis2(pkdtc1))==max(mean(pikis2(pkdtc1)),mean(pikis2(pkdtc2)))
        peaks2=pkdtc1;nonpeaks2=pkdtc2;
    else
        peaks2=pkdtc2;nonpeaks2=pkdtc1;
    end
    idxs=find(peaks2);
    idx=find(nonpeaks);
    peaks(idx(idxs))=true;
    nonpeaks=~peaks;
    cont=cont+1;
   sprintf('peaks:%d',length(find(peaks)))
end

% stores index of the information peaks and its maximum value and the same
% for the noise peaks
newmaxs=maxs(peaks);
newmaxis=maxis(peaks);
noisemaxs=maxs(nonpeaks);
noisemaxis=maxis(nonpeaks);

% in order to eliminate the detected peaks which actually are the same peak
% we compute the median of the difference beteween the peaks (indexes) in
% order to serve as a metric that detects if 2 information peaks are too
% close
mindist=median(diff(newmaxs))/2;
newestmaxs=[];
newestmaxis=[];
%
for(i=1:length(newmaxs))
    %if a certain peak is too close to any other one according to the
    %metric mindist then
   if(length(find(abs(newmaxs-newmaxs(i))<mindist))~=0) 
       % fill the vectors newestmaxs and newestmaxis with the indexes and
       % maximum values of the corrected peaks, i.e. after the near-peak
       % elimination phase
        candidatesmaxis=find(abs(newmaxs-newmaxs(i))<mindist);
        cp=newmaxs(candidatesmaxis);
        candidatesmaxs=newmaxs(candidatesmaxis);
        candidatesmaxis=newmaxis(candidatesmaxis);
        bestcandidate=find(candidatesmaxis==max(candidatesmaxis));
        newestmaxis=[newestmaxis,candidatesmaxis(bestcandidate)];
        newestmaxs=[newestmaxs,candidatesmaxs(bestcandidate)];
   end
end

% given that we repeat the process for the near peaks we have to eliminate
% the repeated corrected peaks *MUDAR ISTO*
finalmaxs=unique(newestmaxs);
finalmaxs
% NEW %
leftmargin=200;
rightmargin=200;
left=finalmaxs-leftmargin;
right=finalmaxs+rightmargin;
if left(1)<1
    left(1)=1;
end
if right(end)>length(sound)
    right(end)=length(sound);
end
energies=zeros(1,length(finalmaxs));
for(i=1:length(energies))
    energies(i)=sum(abs(sound(left(i):right(i)).^2))/(length(left(i):right(i))^2);
end

energies=energies/max(energies);
plot(sound)
for(i=1:length(left))
    line([left(i),left(i)],[-1,1],'color','black')
    line([right(i),right(i)],[-1,1],'color','black')
end







[finalmaxs2,idxs1,idxs2]=unique(newestmaxs);
[logicalidxs,originalidxs]=ismember(finalmaxs,maxs);

finalmaxis=newestmaxis(idxs1);


%initialize cell array that stores the labels of each segment
periodtext=cell(1,length(finalmaxs));


%use the inflection points again in order to distinguish the S1's from the
%S2's
ids=diff(diff(finalmaxis));
difz=diff(sign(ids));
s2=[false,false,difz>0];
s1=[false,false,difz<=0]; 

%fills periodtext accordingly to the inflection points
periodtext(s2)={'S2'};
periodtext(s1)={'S1'};

% special cases: first and last: as the computation of the inflection
% points require 3 consecutive samples in order to compute one value we
% compute the labels of the first and last detected segment
% First
if(finalmaxis(1)>finalmaxis(2))
    periodtext(1)={'S1'};
    periodtext(2)={'S2'};
else
    periodtext(2)={'S1'};
    periodtext(1)={'S2'};
end
%Last
if(finalmaxis(end)>finalmaxis(end-1))
    periodtext(end)={'S1'};
else
     periodtext(end)={'S2'};
end

%computes 2 features: Number of consecutive S1's and the same for the S2's
% the detection of high number of consecutive S1's and S2's suggests the 
% presence of heart pathologies 
consecutiveS1=0;
for(i=1:length(periodtext)-1)
    if(length(setdiff(periodtext(i),{'S1'}))==0 & length(setdiff(periodtext{i+1},{'S1'})))
        consecutiveS1=consecutiveS1+1;
    end
end
consecutiveS2=0;
for(i=1:length(periodtext)-1)
    if(length(setdiff(periodtext(i),{'S2'}))==0 & length(setdiff(periodtext{i+1},{'S2'})))
        consecutiveS2=consecutiveS2+1;
    end
end
consecutives=consecutiveS1+consecutiveS2;
hr=median(diff(finalmaxs));

% plot(sound(1:sound(end)))
% for(i=1:length(finalmaxs))
%     line([finalmaxs(i),finalmaxs(i)],[-1,1],'Color','black')
% end
%Computation of other features like the average of energies of each segment
% the standard deviation etc
avgEnergy=mean(sumis2);
stdEnergy=std(sumis2);
avgDifmaxis=std(finalmaxis);
avgMinis=mean(minis);
avgPikis=mean(pikis);
avgCvars=mean(cvars);
stdCvars=std(cvars);
avgMaxiwavs=mean(maxiwavs);
stdMaxiwavs=std(maxiwavs);
avgdifmaxs=std(diff(finalmaxs));
medianMaxs=median(diff(finalmaxs));

end