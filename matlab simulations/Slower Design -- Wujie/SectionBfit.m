





n=43;
m=352;
SectionBIncreasingField=[zeros(n,1);IncreasingFieldFit(:);zeros(163-n,1)]*30.9;
SectionBCompensationField=[zeros(m,1);CompensationFieldFit(:);zeros(493-m,1)]*115;
SectionBFit=SectionBIncreasingField+SectionBCompensationField;


plot(SectionBFieldFitPosition,SectionBIncreasingField,'black');
plot(SectionBFieldFitPosition,SectionBCompensationField,'green');
plot(SectionBFieldFitPosition,SectionBFit,'blue');


for i=2:1:600; 
    SectionBFieldDeriv(i)=(SectionBFit(i+1)-SectionBFit(i-1))/(SectionBFieldFitPosition(i+1)-SectionBFieldFitPosition(i-1)); % in 1/cm
    AccelaratoinSectionB(i)=(-Bfinal+SectionBFit(i))*SectionBFieldDeriv(i)*(Mu*Lambda)^2;
end

AccelaratoinSectionB(1)=AccelaratoinSectionB(2);
AccelaratoinSectionB(601)=AccelaratoinSectionB(600);
plot(SectionBFieldFitPosition,-AccelaratoinSectionB/amax*100,'red');hold on;

%__________________Theoretical Design



for i=701:1:1225; 
    DesignSectionBFieldDeriv(i)=(DesignFieldProfile(i+1)-DesignFieldProfile(i-1))/(z(i+1)-z(i-1)); % in 1/cm
    DesignAccelaratoinSectionB(i)=(-Bfinal+DesignFieldProfile(i))*DesignSectionBFieldDeriv(i)*(Mu*Lambda)^2;
    SectionBFieldFitPositionB(i)=z(i);
end


plot(SectionBFieldFitPositionB,-DesignAccelaratoinSectionB/amax*100,'red');hold on;


%__________________Practical Design



for i=701:1:1225; 
    AchievedSectionBFieldDeriv(i)=(AchievedFieldProfile(i+1)-AchievedFieldProfile(i-1))/(z(i+1)-z(i-1)); % in 1/cm
    AchievedAccelaratoinSectionB(i)=(-Bfinal+AchievedFieldProfile(i))*AchievedSectionBFieldDeriv(i)*(Mu*Lambda)^2;
    SectionBFieldFitPositionB(i)=z(i);
end


plot(SectionBFieldFitPositionB,-AchievedAccelaratoinSectionB/amax*100,'green');hold on;