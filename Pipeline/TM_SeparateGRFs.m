function [sGRFdata] = TM_SeparateGRFs(MarkerData,GRFdata,Markerset)
%SeparateGRF function separate ground reaction force data into forces onto
%each foot
swindow=10;
warning('on');
%Identifies RTOE, RCAL, LTOE, and LCAL data columns
for i=1:length(Markerset)
    if strcmp(Markerset(i),{'RCAL'})
        rc=i;
    elseif strcmp(Markerset(i),{'RTOE'})
        rt=i;
    elseif strcmp(Markerset(i),{'LCAL'})
        lc=i;
    elseif strcmp(Markerset(i),{'LTOE'})
        lt=i;
    end
end
%finds local maxima of RCAL and LCAL
[~,rclocs,~,rcp]=findpeaks(MarkerData(:,rc*3-1));
for i=2:length(rcp)-1
    if rcp(i)<50
        warning(strcat('low prominence for RCAL peak#_',num2str(i)))
    end
end
[~,lclocs,~,lcp]=findpeaks(MarkerData(:,lc*3-1));
for i=2:length(lcp)-1
    if lcp(i)<50
        warning(strcat('low prominence for LCAL peak#_',num2str(i)))
    end
end
%finds local minima of RTOE and LTOE
[~,rtlocs,~,rtp]=findpeaks(-MarkerData(:,rt*3-1));
for i=2:length(rtp)-1
    if rtp(i)<50
        warning(strcat('low prominence for RTOE peak#_',num2str(i)))
    end
end
[~,ltlocs,~,ltp]=findpeaks(-MarkerData(:,lt*3-1));
for i=2:length(ltp)-1
    if ltp(i)<50
        warning(strcat('low prominence for LTOE peak#_',num2str(i)))
    end
end
%finding nearby gait events according to force plates
p2mratio=(MarkerData(2,1)-MarkerData(1,1))/(GRFdata(2,1)-GRFdata(1,1));
rhsmin=zeros(1,length(rclocs)); rhs=zeros(1,length(rclocs));
for i=1:length(rclocs)
    if rclocs(i)*p2mratio-swindow*p2mratio<=0
        rhsmin(i)=min(GRFdata(1:rclocs(i)*p2mratio+swindow*p2mratio,3));
        rhs(i)=find(GRFdata(1:rclocs(i)*p2mratio+swindow*p2mratio,3)==rhsmin(i),1,'last');
    elseif rclocs(i)*p2mratio+swindow*p2mratio>length(GRFdata)
        rhsmin(i)=min(GRFdata(rclocs(i)*p2mratio-swindow*p2mratio:end,3));
        rhs(i)=find(GRFdata(:,3)==rhsmin(i),1,'last');
    else
        rhsmin(i)=min(GRFdata(rclocs(i)*p2mratio-swindow*p2mratio:rclocs(i)*p2mratio+swindow*p2mratio,3));
        rhs(i)=find(GRFdata(1:rclocs(i)*p2mratio+swindow*p2mratio,3)==rhsmin(i),1,'last');
    end
    if rclocs(i)*p2mratio-swindow*p2mratio>0 && rclocs(i)*p2mratio+swindow*p2mratio<=length(GRFdata)
        if rhsmin(i)>10
            warning(strcat('GRF above 10N (',num2str(rhsmin(i)),'N) at RHS#_',num2str(i)))
        elseif rhsmin(i)==GRFdata(rclocs(i)*p2mratio-swindow*p2mratio,3) || rhsmin(i)==GRFdata(rclocs(i)*p2mratio+swindow*p2mratio,3)
            warning(strcat('minimum GRF on edge of detection window for RHS_#',num2str(i)))
        end
    end
end
lhsmin=zeros(1,length(lclocs)); lhs=zeros(1,length(lclocs));
for i=1:length(lclocs)
    if lclocs(i)*p2mratio-swindow*p2mratio<=0
        lhsmin(i)=min(GRFdata(1:lclocs(i)*p2mratio+swindow*p2mratio,3));
        lhs(i)=find(GRFdata(1:lclocs(i)*p2mratio+swindow*p2mratio,3)==lhsmin(i),1,'last');
    elseif lclocs(i)*p2mratio+swindow*p2mratio>length(GRFdata)
        lhsmin(i)=min(GRFdata(lclocs(i)*p2mratio-swindow*p2mratio:end,3));
        lhs(i)=find(GRFdata(:,3)==lhsmin(i),1,'last');
    else
        lhsmin(i)=min(GRFdata(lclocs(i)*p2mratio-swindow*p2mratio:lclocs(i)*p2mratio+swindow*p2mratio,3));
        lhs(i)=find(GRFdata(1:lclocs(i)*p2mratio+swindow*p2mratio,3)==lhsmin(i),1,'last');
    end
    if lclocs(i)*p2mratio-swindow*p2mratio>0 && lclocs(i)*p2mratio+swindow*p2mratio<=length(GRFdata)
        if lhsmin(i)>10
            warning(strcat('GRF above 10N (',num2str(lhsmin(i)),'N) at LHS#_',num2str(i)))
        elseif lhsmin(i)==GRFdata(lclocs(i)*p2mratio-swindow*p2mratio,3) || lhsmin(i)==GRFdata(lclocs(i)*p2mratio+swindow*p2mratio,3)
            warning(strcat('minimum GRF on edge of detection window for LHS_#',num2str(i)))
        end
    end
end
rtomin=zeros(1,length(rtlocs)); rto=zeros(1,length(rtlocs));
for i=1:length(rtlocs)
    if rtlocs(i)*p2mratio-swindow*p2mratio<=0
        [rtomin(i),rto(i)]=min(GRFdata(1:rtlocs(i)*p2mratio+swindow*p2mratio,12));
    elseif rtlocs(i)*p2mratio+swindow*p2mratio>length(GRFdata)
        [rtomin(i),rto(i)]=min(GRFdata(rtlocs(i)*p2mratio-swindow*p2mratio:end,12));
        rto(i)=rto(i)+rtlocs(i)*p2mratio-swindow*p2mratio-1;
    else
        [rtomin(i),rto(i)]=min(GRFdata(rtlocs(i)*p2mratio-swindow*p2mratio:rtlocs(i)*p2mratio+swindow*p2mratio,12));
        rto(i)=rto(i)+rtlocs(i)*p2mratio-swindow*p2mratio-1;
    end
    if rtlocs(i)*p2mratio-swindow*p2mratio>0 && rtlocs(i)*p2mratio+swindow*p2mratio<=length(GRFdata)
        if rtomin(i)>10
            warning(strcat('GRF above 10N (',num2str(rtomin(i)),'N) at RTO#_',num2str(i)))
        elseif rtomin(i)==GRFdata(rtlocs(i)*p2mratio-swindow*p2mratio,12) || rtomin(i)==GRFdata(rtlocs(i)*p2mratio+swindow*p2mratio,12)
            warning(strcat('minimum GRF on edge of detection window for RTO_#',num2str(i)))
        end
    end
end
ltomin=zeros(1,length(ltlocs)); lto=zeros(1,length(ltlocs));
for i=1:length(ltlocs)
    if ltlocs(i)*p2mratio-swindow*p2mratio<=0
        [ltomin(i),lto(i)]=min(GRFdata(1:ltlocs(i)*p2mratio+swindow*p2mratio,12));
    elseif ltlocs(i)*p2mratio+swindow*p2mratio>length(GRFdata)
        [ltomin(i),lto(i)]=min(GRFdata(ltlocs(i)*p2mratio-swindow*p2mratio:end,12));
        lto(i)=lto(i)+ltlocs(i)*p2mratio-swindow*p2mratio-1;
    else
        [ltomin(i),lto(i)]=min(GRFdata(ltlocs(i)*p2mratio-swindow*p2mratio:ltlocs(i)*p2mratio+swindow*p2mratio,12));
        lto(i)=lto(i)+ltlocs(i)*p2mratio-swindow*p2mratio-1;
    end
    if ltlocs(i)*p2mratio-swindow*p2mratio>0 && ltlocs(i)*p2mratio+swindow*p2mratio<=length(GRFdata)
        if ltomin(i)>10
            warning(strcat('GRF above 10N (',num2str(ltomin(i)),'N) at LTO#_',num2str(i)))
        elseif ltomin(i)==GRFdata(ltlocs(i)*p2mratio-swindow*p2mratio,12) || ltomin(i)==GRFdata(ltlocs(i)*p2mratio+swindow*p2mratio,12)
            warning(strcat('minimum GRF on edge of detection window for LTO_#',num2str(i)))
        end
    end
end

%tests gait events
figure
hold on
plot(MarkerData(1:500,1),MarkerData(1:500,rc*3-1),'k',MarkerData(1:500,1),MarkerData(1:500,lc*3-1),'b')
plot(MarkerData(1:500,1),MarkerData(1:500,rt*3-1),'m',MarkerData(1:500,1),MarkerData(1:500,lt*3-1),'r')
i=1;
while rclocs(i)<500
    plot([MarkerData(rclocs(i),1),MarkerData(rclocs(i),1)],[0,500],'k')
    plot([GRFdata(rhs(i),1),GRFdata(rhs(i),1)],[0,500],'k--')
    i=i+1;
end
i=1;
while lclocs(i)<500
    plot([MarkerData(lclocs(i),1),MarkerData(lclocs(i),1)],[0,500],'b')
    plot([GRFdata(lhs(i),1),GRFdata(lhs(i),1)],[0,500],'b--')
    i=i+1;
end
i=1;
while rtlocs(i)<500
    plot([MarkerData(rtlocs(i),1),MarkerData(rtlocs(i),1)],[0,500],'m')
    plot([GRFdata(rto(i),1),GRFdata(rto(i),1)],[0,500],'m--')
    i=i+1;
end
i=1;
while ltlocs(i)<500
    plot([MarkerData(ltlocs(i),1),MarkerData(ltlocs(i),1)],[0,500],'r')
    plot([GRFdata(lto(i),1),GRFdata(lto(i),1)],[0,500],'r--')
    i=i+1;
end
hold off


%Uses the gait cycle event data to separate GRF
sGRFdata=zeros(size(GRFdata));
sGRFdata(:,1)=GRFdata(:,1);%time data
fweight=GRFdata(:,3)./(GRFdata(:,3)+GRFdata(:,12));
bweight=GRFdata(:,12)./(GRFdata(:,3)+GRFdata(:,12));
%determines starting feature and adds pre-event data
if rhs(1)<lto(1) && rhs(1)<lhs(1) && rhs(1)<rto(1)
    intevent=1;%start rhs, pre-event left single stance
    sGRFdata(1:rhs(1)-1,11:13)=GRFdata(1:rhs(1)-1,2:4)+GRFdata(1:rhs(1)-1,11:13);%combine force for left
    sGRFdata(1:rhs(1)-1,14:16)=GRFdata(1:rhs(1)-1,5:7).*fweight(1:rhs(1)-1)+GRFdata(1:rhs(1)-1,14:16).*bweight(1:rhs(1)-1);%combine COP for left
    sGRFdata(1:rhs(1)-1,17:19)=GRFdata(1:rhs(1)-1,8:10)+GRFdata(1:rhs(1)-1,17:19);%combine moment for left
elseif lto(1)<lhs(1) && lto(1)<rto(1) && lto(1)<rhs(1)
    intevent=2;%start lto, pre-event right led dual stance
    sGRFdata(1:lto(1)-1,2:10)=GRFdata(1:lto(1)-1,2:10);%right on front
    sGRFdata(1:lto(1)-1,11:19)=GRFdata(1:lto(1)-1,11:19);%left on back
elseif lhs(1)<rto(1) && lhs(1)<rhs(1) && lhs(1)<lto(1)
    intevent=3;%start lhs, pre-event right single stance
    sGRFdata(1:lhs(1)-1,2:4)=GRFdata(1:lhs(1)-1,2:4)+GRFdata(1:lhs(1)-1,11:13);%combine force for right
    sGRFdata(1:lhs(1)-1,5:7)=GRFdata(1:lhs(1)-1,5:7).*fweight(1:lhs(1)-1)+GRFdata(1:lhs(1)-1,14:16).*bweight(1:lhs(1)-1);%combine COP for right
    sGRFdata(1:lhs(1)-1,8:10)=GRFdata(1:lhs(1)-1,8:10)+GRFdata(1:lhs(1)-1,17:19);%combine moment for right
elseif rto(1)<rhs(1) && rto(1)<lto(1) && rto(1)<lhs(1)
    intevent=4;%start rto, pre-event left led dual stance
    sGRFdata(1:rto(1)-1,2:10)=GRFdata(1:rto(1)-1,11:19);%right on back
    sGRFdata(1:rto(1)-1,11:19)=GRFdata(1:rto(1)-1,2:10);%left on front
end
%adds data between rhs and lto (right led dual stance)
if intevent==2
    if length(lto)==length(rhs)%starts on lto, ends on rhs
        for i=1:length(rhs)-1
            sGRFdata(rhs(i):lto(i+1)-1,2:10)=GRFdata(rhs(i):lto(i+1)-1,2:10);
            sGRFdata(rhs(i):lto(i+1)-1,11:19)=GRFdata(rhs(i):lto(i+1)-1,11:19);
        end
    else%starts on lto, doesn't end on rhs
        for i=1:length(rhs)
            sGRFdata(rhs(i):lto(i+1)-1,2:10)=GRFdata(rhs(i):lto(i+1)-1,2:10);
            sGRFdata(rhs(i):lto(i+1)-1,11:19)=GRFdata(rhs(i):lto(i+1)-1,11:19);
        end
    end
else
    if length(lto)==length(rhs)%doesn't start on lto, doesn't end on rhs
        for i=1:length(rhs)
            sGRFdata(rhs(i):lto(i)-1,2:10)=GRFdata(rhs(i):lto(i)-1,2:10);
            sGRFdata(rhs(i):lto(i)-1,11:19)=GRFdata(rhs(i):lto(i)-1,11:19);
        end
    else%doesn't start on lto, ends on rhs
        for i=1:length(rhs)-1
            sGRFdata(rhs(i):lto(i)-1,2:10)=GRFdata(rhs(i):lto(i)-1,2:10);
            sGRFdata(rhs(i):lto(i)-1,11:19)=GRFdata(rhs(i):lto(i)-1,11:19);
        end
    end
end
%adds data between lto and lhs (right single stance)
if intevent==3
    if length(lhs)==length(lto)%starts on lhs, ends on lto
        for i=1:length(lto)-1
            sGRFdata(lto(i):lhs(i+1)-1,2:4)=GRFdata(lto(i):lhs(i+1)-1,2:4)+GRFdata(lto(i):lhs(i+1)-1,11:13);
            sGRFdata(lto(i):lhs(i+1)-1,5:7)=GRFdata(lto(i):lhs(i+1)-1,5:7).*fweight(lto(i):lhs(i+1)-1)+GRFdata(lto(i):lhs(i+1)-1,14:16).*bweight(lto(i):lhs(i+1)-1);
            sGRFdata(lto(i):lhs(i+1)-1,8:10)=GRFdata(lto(i):lhs(i+1)-1,8:10)+GRFdata(lto(i):lhs(i+1)-1,17:19);
        end
    else%starts on lhs, doesn't end on lto
        for i=1:length(lto)
            sGRFdata(lto(i):lhs(i+1)-1,2:4)=GRFdata(lto(i):lhs(i+1)-1,2:4)+GRFdata(lto(i):lhs(i+1)-1,11:13);
            sGRFdata(lto(i):lhs(i+1)-1,5:7)=GRFdata(lto(i):lhs(i+1)-1,5:7).*fweight(lto(i):lhs(i+1)-1)+GRFdata(lto(i):lhs(i+1)-1,14:16).*bweight(lto(i):lhs(i+1)-1);
            sGRFdata(lto(i):lhs(i+1)-1,8:10)=GRFdata(lto(i):lhs(i+1)-1,8:10)+GRFdata(lto(i):lhs(i+1)-1,17:19);
        end
    end
else
    if length(lhs)==length(lto)%doesn't start on lhs, doesn't end on lto
        for i=1:length(lto)
            sGRFdata(lto(i):lhs(i)-1,2:4)=GRFdata(lto(i):lhs(i)-1,2:4)+GRFdata(lto(i):lhs(i)-1,11:13);
            sGRFdata(lto(i):lhs(i)-1,5:7)=GRFdata(lto(i):lhs(i)-1,5:7).*fweight(lto(i):lhs(i)-1)+GRFdata(lto(i):lhs(i)-1,14:16).*bweight(lto(i):lhs(i)-1);
            sGRFdata(lto(i):lhs(i)-1,8:10)=GRFdata(lto(i):lhs(i)-1,8:10)+GRFdata(lto(i):lhs(i)-1,17:19);
        end
    else%doesn't start on lhs, ends on lto
        for i=1:length(lto)-1
            sGRFdata(lto(i):lhs(i)-1,2:4)=GRFdata(lto(i):lhs(i)-1,2:4)+GRFdata(lto(i):lhs(i)-1,11:13);
            sGRFdata(lto(i):lhs(i)-1,5:7)=GRFdata(lto(i):lhs(i)-1,5:7).*fweight(lto(i):lhs(i)-1)+GRFdata(lto(i):lhs(i)-1,14:16).*bweight(lto(i):lhs(i)-1);
            sGRFdata(lto(i):lhs(i)-1,8:10)=GRFdata(lto(i):lhs(i)-1,8:10)+GRFdata(lto(i):lhs(i)-1,17:19);
        end
    end
end
%adds data between lhs and rto (left led dual stance)
if intevent==4
    if length(rto)==length(lhs)%starts on rto, ends on lhs
        for i=1:length(lhs)-1
            sGRFdata(lhs(i):rto(i+1)-1,2:10)=GRFdata(lhs(i):rto(i+1)-1,11:19);
            sGRFdata(lhs(i):rto(i+1)-1,11:19)=GRFdata(lhs(i):rto(i+1)-1,2:10);
        end
    else%starts on rto, doesn't end on lhs
        for i=1:length(lhs)
            sGRFdata(lhs(i):rto(i+1)-1,2:10)=GRFdata(lhs(i):rto(i+1)-1,11:19);
            sGRFdata(lhs(i):rto(i+1)-1,11:19)=GRFdata(lhs(i):rto(i+1)-1,2:10);
        end
    end
else
    if length(rto)==length(lhs)%doesn't start on rto, doesn't end on lhs
        for i=1:length(lhs)
            sGRFdata(lhs(i):rto(i)-1,2:10)=GRFdata(lhs(i):rto(i)-1,11:19);
            sGRFdata(lhs(i):rto(i)-1,11:19)=GRFdata(lhs(i):rto(i)-1,2:10);
        end
    else%doesn't start on rto, ends on lhs
        for i=1:length(lhs)-1
            sGRFdata(lhs(i):rto(i)-1,2:10)=GRFdata(lhs(i):rto(i)-1,11:19);
            sGRFdata(lhs(i):rto(i)-1,11:19)=GRFdata(lhs(i):rto(i)-1,2:10);
        end
    end
end
%adds data between rto and rhs (left single stance)
if intevent==1
    if length(rhs)==length(rto)%starts on rhs, ends on rto
        for i=1:length(rto)-1
            sGRFdata(rto(i):rhs(i+1)-1,11:13)=GRFdata(rto(i):rhs(i+1)-1,2:4)+GRFdata(rto(i):rhs(i+1)-1,11:13);
            sGRFdata(rto(i):rhs(i+1)-1,14:16)=GRFdata(rto(i):rhs(i+1)-1,5:7).*fweight(rto(i):rhs(i+1)-1)+GRFdata(rto(i):rhs(i+1)-1,14:16).*bweight(rto(i):rhs(i+1)-1);
            sGRFdata(rto(i):rhs(i+1)-1,17:19)=GRFdata(rto(i):rhs(i+1)-1,8:10)+GRFdata(rto(i):rhs(i+1)-1,17:19);
        end
    else%starts on rhs, doesn't end on rto
        for i=1:length(rto)
            sGRFdata(rto(i):rhs(i+1)-1,11:13)=GRFdata(rto(i):rhs(i+1)-1,2:4)+GRFdata(rto(i):rhs(i+1)-1,11:13);
            sGRFdata(rto(i):rhs(i+1)-1,14:16)=GRFdata(rto(i):rhs(i+1)-1,5:7).*fweight(rto(i):rhs(i+1)-1)+GRFdata(rto(i):rhs(i+1)-1,14:16).*bweight(rto(i):rhs(i+1)-1);
            sGRFdata(rto(i):rhs(i+1)-1,17:19)=GRFdata(rto(i):rhs(i+1)-1,8:10)+GRFdata(rto(i):rhs(i+1)-1,17:19);
        end
    end
else
    if length(rhs)==length(rto)%doesn't start on rhs, doesn't end on rto
        for i=1:length(rto)
            sGRFdata(rto(i):rhs(i)-1,11:13)=GRFdata(rto(i):rhs(i)-1,2:4)+GRFdata(rto(i):rhs(i)-1,11:13);
            sGRFdata(rto(i):rhs(i)-1,14:16)=GRFdata(rto(i):rhs(i)-1,5:7).*fweight(rto(i):rhs(i)-1)+GRFdata(rto(i):rhs(i)-1,14:16).*bweight(rto(i):rhs(i)-1);
            sGRFdata(rto(i):rhs(i)-1,17:19)=GRFdata(rto(i):rhs(i)-1,8:10)+GRFdata(rto(i):rhs(i)-1,17:19);
        end
    else%doesn't start on rhs, ends on rto
        for i=1:length(rto)-1
            sGRFdata(rto(i):rhs(i)-1,11:13)=GRFdata(rto(i):rhs(i)-1,2:4)+GRFdata(rto(i):rhs(i)-1,11:13);
            sGRFdata(rto(i):rhs(i)-1,14:16)=GRFdata(rto(i):rhs(i)-1,5:7).*fweight(rto(i):rhs(i)-1)+GRFdata(rto(i):rhs(i)-1,14:16).*bweight(rto(i):rhs(i)-1);
            sGRFdata(rto(i):rhs(i)-1,17:19)=GRFdata(rto(i):rhs(i)-1,8:10)+GRFdata(rto(i):rhs(i)-1,17:19);
        end
    end
end
%determines last feature and adds post-event data
if rto(end)>rhs(end) && rto(end)>lto(end) && rto(end)>lhs(end)%ends rto, post-event left single stance
    sGRFdata(rto(end):end,11:13)=GRFdata(rto(end):end,2:4)+GRFdata(rto(end):end,11:13);%combine force for left
    sGRFdata(rto(end):end,14:16)=GRFdata(rto(end):end,5:7).*fweight(rto(end):end)+GRFdata(rto(end):end,14:16).*bweight(rto(end):end);%combine COP for left
    sGRFdata(rto(end):end,17:19)=GRFdata(rto(end):end,8:10)+GRFdata(rto(end):end,17:19);%combine moment for left
elseif rhs(end)>lto(end) && rhs(end)>lhs(end) && rhs(end)>rto(end)%ends rhs, post-event right led dual stance
    sGRFdata(rhs(end):end,2:10)=GRFdata(rhs(end):end,2:10);%right on front
    sGRFdata(rhs(end):end,11:19)=GRFdata(rhs(end):end,11:19);%left on back
elseif lto(end)>lhs(end) && lto(end)>rto(end) && lto(end)>rhs(end)%ends lto, post-event right single stance
    sGRFdata(lto(end):end,2:4)=GRFdata(lto(end):end,2:4)+GRFdata(lto(end):end,11:13);%combine force for right
    sGRFdata(lto(end):end,5:7)=GRFdata(lto(end):end,5:7).*fweight(lto(end):end)+GRFdata(lto(end):end,14:16).*bweight(lto(end):end);%combine COP for right
    sGRFdata(lto(end):end,8:10)=GRFdata(lto(end):end,8:10)+GRFdata(lto(end):end,17:19);%combine moment for right
elseif lhs(end)>rto(end) && lhs(end)>rhs(end) && lhs(end)>lto(end)%ends lhs, post-event left led dual stance
    sGRFdata(lhs(end):end,2:10)=GRFdata(lhs(end):end,11:19);%right on back
    sGRFdata(lhs(end):end,11:19)=GRFdata(lhs(end):end,2:10);%left on front
end
for i=1:floor(length(MarkerData)/2000)
    figure
    hold on
    subplot(2,1,1)
    plot(GRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,1),GRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,3),GRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,1),GRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,12))
    subplot(2,1,2)
    plot(sGRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,1),sGRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,3),sGRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,1),sGRFdata((i-1)*2000*p2mratio+1:i*2000*p2mratio,12))
end
figure
hold on
subplot(2,1,1)
plot(GRFdata(i*2000*p2mratio+1:end,1),GRFdata(i*2000*p2mratio+1:end,3),GRFdata(i*2000*p2mratio+1:end,1),GRFdata(i*2000*p2mratio+1:end,12))
subplot(2,1,2)
plot(sGRFdata(i*2000*p2mratio+1:end,1),sGRFdata(i*2000*p2mratio+1:end,3),sGRFdata(i*2000*p2mratio+1:end,1),sGRFdata(i*2000*p2mratio+1:end,12))
end