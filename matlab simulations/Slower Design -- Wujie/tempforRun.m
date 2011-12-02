shift=50;


for i=701:1:1225; 
    AchievedSectionBFieldDeriv(i)=(AchievedFieldProfile(i+1)-AchievedFieldProfile(i-1))/(z(i+1)-z(i-1)); % in 1/cm
    AchievedAccelaratoinSectionB(i)=(-Bfinal+AchievedFieldProfile(i))*AchievedSectionBFieldDeriv(i)*(Mu*Lambda)^2;
    SectionBFieldFitPositionB(i)=z(i);
end

figure(2);
plot(SectionBFieldFitPositionB,-AchievedAccelaratoinSectionB/amax*100,'green');hold on;


shift=50;

for i=1:1:562; 
    AchievedDecreasingFieldDeriv(i)=(AchievedFieldProfile(i+1+shift)-AchievedFieldProfile(i-1+shift))/(z(i+1+shift)-z(i-1+shift)); % in 1/cm
    AchievedAccelaratoinSectionA(i)=(-Bfinal+AchievedFieldProfile(i+shift))*AchievedDecreasingFieldDeriv(i)*(Mu*Lambda)^2;
    DecreasingFieldFitPosition(i)=z(i+50);
end


plot(DecreasingFieldFitPosition,-AchievedAccelaratoinSectionA/amax*100,'green');