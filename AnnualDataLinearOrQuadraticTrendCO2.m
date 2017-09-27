%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads annual mean CO2 data from Cape Grim, then does OLS and   %
% GLS for both linear (y = b0 + b1.t) and quadratic (y = b0 + b1.t + b2.t^2). %
% This code is for Matlab or Octave (no special libraries required).          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doQuatratic=false;  % otherwise linear fit

plotACF=false;      % plot autocorrelation function of residuals?

plotLongTerm=4;     % plot long-term extrapolated model with 95% conf. int.: 0=no, 1=OLS, 2=ACF GLS, 3=DW GLS, 4=new GLS

disp('========================================================');

% Read data as csv numeric file of two columns: year (AD) and CO2 (ppm)
% This file is from data at http://www.csiro.au/greenhouse-gases/
CO2data=csvread('CapeGrimAnnualCO2.csv');
disp('Data: Atmospheric background carbon dioxide at Cape Grim');
disp('from http://www.csiro.au/greenhouse-gases/');

%%%%%%%%%%%%%%%%%%%%%
% Assemble the data %
%%%%%%%%%%%%%%%%%%%%%
[n,~]=size(CO2data);
t=0;
x=0;
y=0;
for i=1:n
    t(i)=CO2data(i,1);
    x(i,1)=1;
    x(i,2)=t(i)-2000;
    y(i,1)=CO2data(i,2);
    if doQuatratic
        x(i,3)=x(i,2)*x(i,2);
    end
end

disp('--------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,k]=size(x);
DegreesOfFreedom=n-k;
varName(1,:)='(Intercept)';
varName(2,:)='year - 2000';
if doQuatratic
    varName(3,:)='(yr-2000)^2';
end
disp('OLS calculation by matrix algebra:');
betaOLS=inv(transpose(x)*x)*transpose(x)*y;
residuals=y-x*betaOLS;

SSR=transpose(residuals)*residuals;
% for improved rho estimate:
meanPred=(1/n)*ones(n,n)*x*betaOLS; % n x 1 vector all elements=mean predicted value
predDev=x*betaOLS-meanPred;
SSPD=transpose(predDev)*predDev;

residStdErrOLS=std(residuals)*sqrt((n-1)/(n-k));
stdErrorOLS=0;
tValue=0;
normCov=inv(transpose(x)*x);
for m=1:k
    stdErrorOLS(m,1)=residStdErrOLS*sqrt(normCov(m,m));
    tValue(m,1)=betaOLS(m,1)/stdErrorOLS(m,1);
end
disp('                 Estimate     Std. Error    t value');
for m=1:k
    fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaOLS(m,1),stdErrorOLS(m,1),tValue(m,1));
end
disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k),' degrees of freedom']);
rmsError=sqrt((transpose(residuals)*residuals)/n);
disp(['Model RMS error: ',num2str(rmsError),' (ppm CO2)']);

if doQuatratic
    for i=1:n
        bestFitTunEq(i)=t(i);
        bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1)+x(i,3)*betaOLS(3,1);
    end
else
    for i=1:n
        bestFitTunEq(i)=t(i);
        bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1);
    end
end


disp('--------------------------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and plot sample autocorrelation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample autocorrelation function (provided in MatLab but not Octave)
Ymean=mean(residuals);
R=0;
Lag=0;
zeroLine=0;
CIu=0;
CIl=0;
denom=0;
for i=1:n
    denom=denom+(residuals(i)-Ymean)^2;
end
for ka=0:20
    numer=0;
    for na=1:n-ka
        numer=numer+(residuals(na)-Ymean)*(residuals(na+ka)-Ymean);
    end
    R(ka+1,1)=numer/denom;
    Lag(ka+1)=ka;
    zeroLine(ka+1)=0;
    CIu(ka+1)=1.96/sqrt(n);
    CIl(ka+1)=-1.96/sqrt(n);
end
%disp(['Residual ACF lag 1 = ',num2str(R(2)),' (first estimate of rho)']);
if plotACF
    plot(Lag,R,'.-b',Lag,zeroLine,'-g',Lag,CIu,'-r',Lag,CIl,'-r');
    axis([0 20 -0.4 1]);
    xlabel('Lag');
    ylabel('Sample Autocorrelation');
    grid('on');
    if doQuatratic
        title(['Cape Grim CO_2 annual data, ACF of quadratic model residuals']);
    else
        title(['Cape Grim CO_2 annual data, ACF of linear model residuals']);
    end
    pause;
end

rho=R(2,1); % for GLS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple n_eff = n(1-rho)/(1+rho) method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t = b*sqrt(n-k)/sqrt(SSR/sum of squares of x differences from mean)
ne=n*(1-rho)/(1+rho);
tScale=sqrt((1-rho)/(1+rho));
disp(['n_eff = n(1-rho)/(1+rho) = ',num2str(ne),' -> t-slope = ',num2str(tValue(2,1)*sqrt((ne-k)/(n-k)))]);
disp(['or scale OLS t values by sqrt[(1-rho)/(1+rho)] = ',num2str(tScale)]);
disp(['                 -> t-slope = ',num2str(tValue(2,1)*tScale)]);
if doQuatratic
    disp(['                 -> t-curve = ',num2str(tValue(3,1)*tScale)]);
end
disp('--------------------------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['GLS using ACF lag 1 for AR1 rho = ',num2str(rho)]);

% Conventional GLS for AR1 process, but calculate Sinv directly:
% initialise Sinv with diagonal elements = (1+rho^2)/(1-rho^2)
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
% change the two corner elements to 1/(1-rho^2)
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
% now fill in the super and sub diagonal elements as -rho/(1-rho^2)
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
% Note that with the above formulation of Sinv, all diagonal elements of inv(Sinv) = 1.

betaGLSacf=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSacf;
residStdErrACF=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
stdErrorACF=0;
tValue=0;
normCov=inv(transpose(x)*Sinv*x);
for m=1:k
    stdErrorACF(m,1)=residStdErrACF*sqrt(normCov(m,m));
    tValue(m,1)=betaGLSacf(m,1)/stdErrorACF(m,1);
end
disp('t-value = beta / std.error');
disp('                 Estimate     Std. Error    t value');
for m=1:k
    fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSacf(m,1),stdErrorACF(m,1),tValue(m,1));
end
disp(['Residual standard error: ',num2str(residStdErrACF),' on ',num2str(n-k),' degrees of freedom']);
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Model RMS error: ',num2str(rmsError),' (ppm CO2)']);

if doQuatratic
    for i=1:n
        bestFitTacf(i)=t(i);
        bestFitYacf(i)=betaGLSacf(1,1)+x(i,2)*betaGLSacf(2,1)+x(i,3)*betaGLSacf(3,1);
    end
else
    for i=1:n
        bestFitTacf(i)=t(i);
        bestFitYacf(i)=betaGLSacf(1,1)+x(i,2)*betaGLSacf(2,1);
    end
end



disp('--------------------------------------------------------');

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Durbin-Watson Statistic, rho estimate and improved estimate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSRFD=0;
for i=2:n
    SSRFD=SSRFD+(residuals(i)-residuals(i-1))^2;
end
dw=SSRFD/SSR;

% rho from Durbin-Watson
rhoDW=1-0.5*dw;

% rho from modified Durbin-Watson
rhoDWM=1-0.5*SSRFD/(SSR*(2-exp(-SSPD/SSR)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=rhoDW;
disp(['GLS using Durbin-Watson for AR1 rho = ',num2str(rho)]);
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
betaGLSdw=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSdw;
residStdErrDW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
stdErrorDW=0;
tValue=0;
normCov=inv(transpose(x)*Sinv*x);
for m=1:k
    stdErrorDW(m,1)=residStdErrDW*sqrt(normCov(m,m));
    tValue(m,1)=betaGLSdw(m,1)/stdErrorDW(m,1);
end
disp('t-value = beta / std.error');
disp('                 Estimate     Std. Error    t value');
for m=1:k
    fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSdw(m,1),stdErrorDW(m,1),tValue(m,1));
end
disp(['Residual standard error: ',num2str(residStdErrDW),' on ',num2str(n-k),' degrees of freedom']);
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Model RMS error: ',num2str(rmsError),' (ppm CO2)']);

if doQuatratic
    for i=1:n
        bestFitTdw(i)=t(i);
        bestFitYdw(i)=betaGLSdw(1,1)+x(i,2)*betaGLSdw(2,1)+x(i,3)*betaGLSdw(3,1);
    end
else
    for i=1:n
        bestFitTdw(i)=t(i);
        bestFitYdw(i)=betaGLSdw(1,1)+x(i,2)*betaGLSdw(2,1);
    end
end


disp('--------------------------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=rhoDWM;
disp(['GLS using modified Durbin-Watson for AR1 rho = ',num2str(rho)]);
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
betaGLSmdw=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSmdw;
residStdErrMDW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
stdErrorMDW=0;
tValue=0;
normCov=inv(transpose(x)*Sinv*x);
for m=1:k
    stdErrorMDW(m,1)=residStdErrMDW*sqrt(normCov(m,m));
    tValue(m,1)=betaGLSmdw(m,1)/stdErrorMDW(m,1);
end
disp('t-value = beta / std.error');
disp('                 Estimate     Std. Error    t value');
for m=1:k
    fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSmdw(m,1),stdErrorMDW(m,1),tValue(m,1));
end
disp(['Residual standard error: ',num2str(residStdErrMDW),' on ',num2str(n-k),' degrees of freedom']);
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Model RMS error: ',num2str(rmsError),' (ppm CO2)']);

if doQuatratic
    for i=1:n
        bestFitTmdw(i)=t(i);
        bestFitYmdw(i)=betaGLSmdw(1,1)+x(i,2)*betaGLSmdw(2,1)+x(i,3)*betaGLSmdw(3,1);
    end
else
     for i=1:n
        bestFitTmdw(i)=t(i);
        bestFitYmdw(i)=betaGLSmdw(1,1)+x(i,2)*betaGLSmdw(2,1);
    end
end


disp('========================================================');


% Plot regression lines:
plot(t,y,'k.-',bestFitTunEq,bestFitYunEq,'r*:',bestFitTacf,bestFitYacf,'gx-',bestFitTdw,bestFitYdw,'bo-',bestFitTmdw,bestFitYmdw,'m+-');
legend('observed','OLS fit',['GLS ACF lag-1 AR1(',num2str(R(2),3),')'],['GLS Durbin-Watson AR1(',num2str(rhoDW,3),')'],['GLS Modified D-W AR1(',num2str(rhoDWM,3),')'],'Location','northwest');
title([num2str(n),' measured data points, and best fit lines, Cape Grim CO_2 annual data']);
xlabel('Year');
ylabel('CO_2 Concentration, ppm');
if plotLongTerm>0.5
    pause;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output long-term projection 1960 - 2040 as required %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include uncertainty of the regression coefficients as well as the residual standard error
% to estimate confidence interval beyond the measured data.
if doQuatratic
    for i=1:81
        yr=i-41;
        yx=yr+2000-1997;  % years away from the centre of the data
        YearLT(i)=yr+2000; % year AD
        OLSLT(i)=betaOLS(1,1)    +betaOLS(2,1)*yr    +betaOLS(3,1)*yr*yr;
        OLSLTci(i)=1.96*sqrt(residStdErrOLS^2+(yx*stdErrorOLS(2,1))^2+(yx*yx*stdErrorOLS(3,1))^2); % +/- 95% CI
        ACFLT(i)=betaGLSacf(1,1) +betaGLSacf(2,1)*yr +betaGLSacf(3,1)*yr*yr;
        ACFLTci(i)=1.96*sqrt(residStdErrACF^2+(yx*stdErrorACF(2,1))^2+(yx*yx*stdErrorACF(3,1))^2); % +/- 95% CI
        GLSdwLT(i)=betaGLSdw(1,1)  +betaGLSdw(2,1)*yr  +betaGLSdw(3,1)*yr*yr;
        GLSdwLTci(i)=1.96*sqrt(residStdErrDW^2+(yx*stdErrorDW(2,1))^2+(yx*yx*stdErrorDW(3,1))^2); % +/- 95% CI
        GLSmdwLT(i)=betaGLSmdw(1,1)+betaGLSmdw(2,1)*yr+betaGLSmdw(3,1)*yr*yr;
        GLSmdwLTci(i)=1.96*sqrt(residStdErrMDW^2+(yx*stdErrorMDW(2,1))^2+(yx*yx*stdErrorMDW(3,1))^2); % +/- 95% CI
    end

    if plotLongTerm==1
        plot(t,y,'k-',YearLT,OLSLT+OLSLTci,'m-.',YearLT,OLSLT,'r-',YearLT,OLSLT-OLSLTci,'b-.');
        legend('observed','95% conf. int. upper','quadratic OLS fit','model 95% confidence interval lower','Location','northwest');
    elseif plotLongTerm==2
        plot(t,y,'k-',YearLT,ACFLT+ACFLTci,'m-.',YearLT,ACFLT,'r-',YearLT,ACFLT-ACFLTci,'b-.');
        legend('observed','model 95% confidence interval upper','quadratic GLS ACF lag-1 model','model 95% confidence interval lower','Location','northwest');
    elseif plotLongTerm==3
        plot(t,y,'k-',YearLT,GLSdwLT+GLSdwLTci,'m-.',YearLT,GLSdwLT,'r-',YearLT,GLSdwLT-GLSdwLTci,'b-.');
        legend('observed','model 95% confidence interval upper','quadratic GLS Durbin-Watson','model 95% confidence interval lower','Location','northwest'); 
    elseif plotLongTerm==4
        plot(t,y,'k-',YearLT,GLSmdwLT+GLSmdwLTci,'m-.',YearLT,GLSmdwLT,'r-',YearLT,GLSmdwLT-GLSmdwLTci,'b-.');
        legend('observed','model 95% confidence interval upper','quadratic GLS modified Durbin-Watson','model 95% confidence interval lower','Location','northwest');
    end
    if plotLongTerm>0.5
        xlabel('Year');
        ylabel('CO_2 Concentration, ppm');
        axis([1960 2040 300 440]);
        title([num2str(n),' measured data points, and fitted model, extrapolated, with confidence interval, Cape Grim CO_2']);
    end
else
    for i=1:81
        yr=i-41;
        yx=yr+2000-1997;  % years away from the centre of the data
        YearLT(i)=yr+2000; % year AD
        OLSLT(i)=betaOLS(1,1)    +betaOLS(2,1)*yr;
        OLSLTci(i)=1.96*sqrt(residStdErrOLS^2+(yx*stdErrorOLS(2,1))^2); % +/- 95% CI
        ACFLT(i)=betaGLSacf(1,1) +betaGLSacf(2,1)*yr;
        ACFLTci(i)=1.96*sqrt(residStdErrACF^2+(yx*stdErrorACF(2,1))^2); % +/- 95% CI
        GLSdwLT(i)=betaGLSdw(1,1)  +betaGLSdw(2,1)*yr;
        GLSdwLTci(i)=1.96*sqrt(residStdErrDW^2+(yx*stdErrorDW(2,1))^2); % +/- 95% CI
        GLSmdwLT(i)=betaGLSmdw(1,1)+betaGLSmdw(2,1)*yr;
        GLSmdwLTci(i)=1.96*sqrt(residStdErrMDW^2+(yx*stdErrorMDW(2,1))^2); % +/- 95% CI
    end

    if plotLongTerm==1
        plot(t,y,'k-',YearLT,OLSLT+OLSLTci,'m-.',YearLT,OLSLT,'r-',YearLT,OLSLT-OLSLTci,'b-.');
        legend('observed','model 95% confidence interval upper','linear OLS fit','model 95% confidence interval lower','Location','northwest');
    elseif plotLongTerm==2
        plot(t,y,'k-',YearLT,ACFLT+ACFLTci,'m-.',YearLT,ACFLT,'r-',YearLT,ACFLT-ACFLTci,'b-.');
        legend('observed','model 95% confidence interval upper','linear GLS ACF lag-1 model','model 95% confidence interval lower','Location','northwest');
    elseif plotLongTerm==3
        plot(t,y,'k-',YearLT,GLSdwLT+GLSdwLTci,'m-.',YearLT,GLSdwLT,'r-',YearLT,GLSdwLT-GLSdwLTci,'b-.');
        legend('observed','model 95% confidence interval upper','linear GLS Durbin-Watson','model 95% confidence interval lower','Location','northwest'); 
    elseif plotLongTerm==4
        plot(t,y,'k-',YearLT,GLSmdwLT+GLSmdwLTci,'m-.',YearLT,GLSmdwLT,'r-',YearLT,GLSmdwLT-GLSmdwLTci,'b-.');
        legend('observed','model 95% confidence interval upper','linear GLS modified Durbin-Watson','model 95% confidence interval lower','Location','northwest');
    end
    if plotLongTerm>0.5
        xlabel('Year');
        ylabel('CO_2 Concentration, ppm');
        axis([1960 2040 300 440]);
        title([num2str(n),' measured data points, and fitted model, extrapolated, with confidence interval, Cape Grim CO_2']);
    end
end

