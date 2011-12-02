%__________________Theoretical Design

shift=100;

for i=1:1:562; 
    DesignDecreasingFieldDeriv(i)=(DesignFieldProfile(i+1+shift)-DesignFieldProfile(i-1+shift))/(z(i+1+shift)-z(i-1+shift)); % in 1/cm
    DesignAccelaratoinSectionA(i)=(-Bfinal+DesignFieldProfile(i+shift))*DesignDecreasingFieldDeriv(i)*(Mu*Lambda)^2;
    DecreasingFieldFitPosition(i)=z(i+shift);
end


plot(DecreasingFieldFitPosition,-DesignAccelaratoinSectionA/amax*100,'red');hold on;


%__________________Practical Design

shift=100;

for i=1:1:562; 
    AchievedDecreasingFieldDeriv(i)=(AchievedFieldProfile(i+1+shift)-AchievedFieldProfile(i-1+shift))/(z(i+1+shift)-z(i-1+shift)); % in 1/cm
    AchievedAccelaratoinSectionA(i)=(-Bfinal+AchievedFieldProfile(i+shift))*AchievedDecreasingFieldDeriv(i)*(Mu*Lambda)^2;
    DecreasingFieldFitPosition(i)=z(i+shift);
end


plot(DecreasingFieldFitPosition,-AchievedAccelaratoinSectionA/amax*100,'green');


%_________________Practical Achieved




for i=2:1:561; 
    DecreasingFieldDeriv(i)=(DecreasingFieldFit(i+1)-DecreasingFieldFit(i-1))/(DecreasingFieldFitPosition(i+1)-DecreasingFieldFitPosition(i-1)); % in 1/cm
    AccelaratoinSectionA(i)=(-Bfinal+DecreasingFieldFit(i))*DecreasingFieldDeriv(i)*(Mu*Lambda)^2;
end

DecreasingFieldDeriv(562)=(DecreasingFieldFit(562)-DecreasingFieldFit(561))/0.1;



AccelaratoinSectionA(562)=(-Bfinal+DecreasingFieldFit(562))*DecreasingFieldDeriv(562)*(Mu*Lambda)^2;

AccelaratoinSectionA(1)=0;

plot(DecreasingFieldFitPosition,-AccelaratoinSectionA/amax*100,'blue');
