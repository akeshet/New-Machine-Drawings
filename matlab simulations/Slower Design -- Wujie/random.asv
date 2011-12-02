
tempField=FitMField1;

temppositionB=positionB-40;
tempField=tempField+bfield1(FitMPosition, -1*3, tube_OD, temppositionB-TurnsB(1)*wire_thickness-wire_thickness*(2*WeakTurnsB(1)+3), 3, -wire_thickness);

plot(FitMPosition,tempField)


figure(4)

n=980;
tempField=FitFinalFieldTotal-[zeros(n,1);FitComMField(1:1041-n)*75];


plot(FitFinalFieldPosition,tempField)




figure(4)

tempField=StimulationField+bfield1(StimulationPosition, 1*30, tube_OD+16*wire_thickness, wire_thickness*(TurnsA(10)+1), 2, wire_thickness);

plot(StimulationPosition,tempField)


n=298

temmp=[FitMField(:);zeros(20,1)]+[zeros(n,1);-115*FitComMField(:);zeros(312-n,1)];
plot(FitMPosition,temmp);



qq=55:0.1:61;
BeforeMPosition=[DecreasingFieldFitPosition(:);qq(:);FitMPosition(:)];
BeforeMField=[DecreasingFieldFit(:);zeros(61,1);temmp(:)];
plot(BeforeMPosition,BeforeMField);


BeforeFinalBField=[DecreasingFieldFit(:);zeros(10,1);]