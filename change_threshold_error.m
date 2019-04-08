% defining a small (threshold) error
results_new=0.0;
%xlswrite('Results_under_threshold.xls',results_new)
threshold_error=getinfo7;
% finding solutions for 'good fits' (under threshold) for the thermograms
results_new=results(results(:,nspec+3)<=threshold_error,:);
xlswrite('Results_under_threshold.xls',results_new)
% defining the weighting factor
inverse_error=1./results_new(:,nspec+3);
% calculating the averaged estimates of the properties
for qd=1:nspec
average_X(qd)=(sum(results_new(:,qd).'*inverse_error(:)))/sum(inverse_error);
end
% average_X1=(sum(results_new(:,1).'*inverse_error(:)))/sum(inverse_error);
% average_X2=(sum(results_new(:,2).'*inverse_error(:)))/sum(inverse_error);
% average_X3=(sum(results_new(:,3).'*inverse_error(:)))/sum(inverse_error);
% average_X4=(sum(results_new(:,4).'*inverse_error(:)))/sum(inverse_error);
average_dHvap=(sum(results_new(:,nspec+1).'*inverse_error(:)))/sum(inverse_error);
average_alpha=(sum(log10(results_new(:,nspec+2)).'*inverse_error(:)))/sum(inverse_error);
% calculating the averaged estimated thermogram and the standard deviation
% (uncertainty) of the estimated thermogram
for qe=1:ntrials
average_MFR(qe)=(sum(results_new(:,nspec+3+qe).'*inverse_error(:)))/sum(inverse_error);
%stdev_MFR(qe)=sqrt(sum((results_new(:,nspec+3+qe).'-average_MFR(qe)).^2.*inverse_error(:)')/sum(inverse_error));
%stdev_MFR(qe)=0.001;
stdev_MFR_max(qe)=max(results_new(:,nspec+3+qe));
stdev_MFR_min(qe)=min(results_new(:,nspec+3+qe));
%stdev_MFR_max(qe)=average_MFR(qe);
%stdev_MFR_min(qe)=average_MFR(qe);
end
% calculating the standard deviation(uncertainty) of the properties
for qf=1:nspec
stdev_X(qf)=sqrt(sum((results_new(:,qf).'-average_X(qf)).^2.*inverse_error(:)')/sum(inverse_error));
end
% stdev_X1=sqrt(sum((results_new(:,1).'-average_X1).^2.*inverse_error(:)')/sum(inverse_error));
% stdev_X2=sqrt(sum((results_new(:,2).'-average_X2).^2.*inverse_error(:)')/sum(inverse_error));
% stdev_X3=sqrt(sum((results_new(:,3).'-average_X3).^2.*inverse_error(:)')/sum(inverse_error));
% stdev_X4=sqrt(sum((results_new(:,4).'-average_X4).^2.*inverse_error(:)')/sum(inverse_error));
stdev_dHvap=sqrt(sum((results_new(:,nspec+1).'-average_dHvap).^2.*inverse_error(:)')/sum(inverse_error));
stdev_alpha=sqrt(sum((log10(results_new(:,nspec+2)).'-average_alpha).^2.*inverse_error(:)')/sum(inverse_error));
%%%%%%%%%%%%%%%% Making%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Volatility distribution
figure(1)
xbar=(1:nspec);
ybar=(average_X);
bar(xbar,ybar,'r')
xlabel('C^* (\mug m^-^3)')
%set(gca,'xtick',1:length(ybar))
set(gca,'Xtick',1:nspec,'xticklabel',cstar)
ylabel('Mass Fraction')
hold on
errorbar(xbar, ybar, (stdev_X),'+k')
set(gca,'yLim', [0 1])
sdir=['./','Volatility Distribution','.','jpg'];
saveas(gcf,sdir,'tif');
saveas(gcf,sdir(1:end-4),'fig')
%print(gcf,'-djpeg','-r1000','Volatility Distribution')
%%%
% Vaporization enthalpy
figure(2)
xbar=[4 5 6];
ybar=[0 average_dHvap./1000. 0];
bar(xbar,ybar,'m')
ylabel('\DeltaH_v_a_p (kJ mol ^-^1)')
set(gca,'xtick',1:length(ybar))
set(gca,'xticklabel','')
ylim([0 140])
hold on
errorbar(xbar,ybar,[0 stdev_dHvap./1000. 0],'+k')
sdir=['./','dHvap','.','jpg'];
saveas(gcf,sdir,'tif');
saveas(gcf,sdir(1:end-4),'fig')
% Accommodation coefficient
figure(3)
xbar=[-2 -1 0];
ybar=[0 10.^(average_alpha) 0];
L=[0 abs(10.^(average_alpha)-10.^(average_alpha-stdev_alpha)) 0];
U=[0 abs(10.^(average_alpha)-10.^(average_alpha+stdev_alpha)) 0];
bar(xbar,ybar,'c', 'baseValue',0.00001)
ylabel('Accommodation coefficient')
set(gca,'xtick',1:length(ybar))
set(gca,'xticklabel','')
ylim([0.00001 2])
hold on
errorbar(xbar,ybar,L,U,'+k') 
set(gca,'YScale','log')
sdir=['./','Accomodation coefficient','.','jpg'];
saveas(gcf,sdir,'tif');
saveas(gcf,sdir(1:end-4),'fig')
% Thermograms
% Area plot for thermograms from experiments (with red dots) and the
% predicted with an area of uncertainty +/-1 sigma
figure(4)
x=T_f-273.15; 
y=average_MFR;
%e=stdev_MFR; % error=stdev
e1=stdev_MFR_max-average_MFR;
e2=average_MFR-stdev_MFR_min;
x(isinf(y)==1)=[];
%e(isinf(y)==1)=[];
e1(isinf(y)==1)=[];
e2(isinf(y)==1)=[];
y(isinf(y)==1)=[];
%Y=[y'-e',2*e']; % [value-stdev,2*stdev]
Y=[y'-e2',e2'+e1'];
h=area(x,Y);
hold on
plot(x,y,'-k','LineWidth',3)
hold off
set(h(1),'FaceColor',[1 1 1])
set(h(1),'EdgeColor',[1 1 1])
set(h(2),'FaceColor',[0.9 0.9 0.9])
set(gca,'TickDir','out','XMinortick','on','YMinortick','on',...
         'Layer','top') % Properties of the axes, gca: get current axes
set(gca,'XLim',[20 x(end)])
set(gca,'XTick',[20:20:x(end)])
%set(gca,'YLim',[0.0 1.0])
hold on
plot(x,experimental,'or','MarkerSize',10,'MarkerFaceColor',[1 0 0])
xlabel('Temperature (^oC)', 'LineWidth', 14)
ylabel('MFR', 'LineWidth', 14)
title('Thermogram')
sdir=['./','Thermograms','.','jpg'];
saveas(gcf,sdir,'tif');
saveas(gcf,sdir(1:end-4),'fig')
% Save results of estimation of averages and uncertainties (standard
% deviations)in excel file. Specifically for volatility distribution, 
%vaporization enthalpy, accommodation coefficient and predicted thermograms 
%(MFRs). The excel file has the name, when saved, 'Estimations and
%Uncertainties'
row_header1={'C*(ug/m3)'};
row_header2={'Avg_X'};
row_header3={'Stdev_X'};
row_header4={'Avg_dHvap(kJ/mol)'};
row_header5={'Stdev_dHvap'};
row_header6={'Avg_am'};
row_header7={'+sigma'};
row_header8={'-sigma'};
row_header9={'Avg_MFR'};
%row_header10={'Stdev_MFR'};
row_header10={'Max_MFR'};
row_header11={'Min_MFR'};
xlswrite('Estimations and Uncertainties',row_header1,'Sheet1','A2');
xlswrite('Estimations and Uncertainties',row_header2,'Sheet1','A3');
xlswrite('Estimations and Uncertainties',row_header3,'Sheet1','A4');
xlswrite('Estimations and Uncertainties',row_header4,'Sheet1','A6');
xlswrite('Estimations and Uncertainties',row_header5,'Sheet1','A7');
xlswrite('Estimations and Uncertainties',row_header6,'Sheet1','A9');
xlswrite('Estimations and Uncertainties',row_header7,'Sheet1','A10');
xlswrite('Estimations and Uncertainties',row_header8,'Sheet1','A11');
xlswrite('Estimations and Uncertainties',row_header9,'Sheet1','A13');
xlswrite('Estimations and Uncertainties',row_header10,'Sheet1','A14');
xlswrite('Estimations and Uncertainties',row_header11,'Sheet1','A15');
xlswrite('Estimations and Uncertainties',cstar,'Sheet1','B2')
xlswrite('Estimations and Uncertainties',average_X,'Sheet1','B3')
xlswrite('Estimations and Uncertainties',stdev_X,'Sheet1','B4')
xlswrite('Estimations and Uncertainties',average_dHvap/1000.,'Sheet1','B6')
xlswrite('Estimations and Uncertainties',stdev_dHvap/1000.,'Sheet1','B7')
xlswrite('Estimations and Uncertainties',10.^(average_alpha),'Sheet1','B9')
xlswrite('Estimations and Uncertainties',U(2),'Sheet1','B10')
xlswrite('Estimations and Uncertainties',L(2),'Sheet1','B11')
xlswrite('Estimations and Uncertainties',average_MFR,'Sheet1','B13')
%xlswrite('Estimations and Uncertainties',stdev_MFR,'Sheet1','B14')
xlswrite('Estimations and Uncertainties',stdev_MFR_max,'Sheet1','B14')
xlswrite('Estimations and Uncertainties',stdev_MFR_min,'Sheet1','B15')
%save all results from all combinations of the properties that were tested
xlswrite('Results.xls', results)
% save all errors (fits) from the combinations that were tested
xlswrite('Results_error.xls',error')