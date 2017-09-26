%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads annual temperature anomaly data,               %
% then does OLS and GLS for a linear (y = b0 + b1.t) model.         %
% This code is for Matlab or Octave (no special libraries required).%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotACF=false;      % plot autocorrelation function of residuals?

disp('========================================================');

% Read data as csv numeric file of two columns: year (AD) and reconstructed temperature anomaly (K)
% This file is from data at https://www1.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/antarctica-temp2006.txt
% but only including the data from 1850 to 1999.
deltaTdata=csvread('AntarcticaAnnualDT.csv');
disp('Data: Antarctic temperature anomaly (last 150 years analysed).');
disp('Schneider, Steig, van Ommen, Dixon, Mayewski, Jones, and Bitz 2006');



%%%%%%%%%%%%%%%%%%%%%
% Assemble the data %
%%%%%%%%%%%%%%%%%%%%%
[n,~]=size(deltaTdata);
t=0;
x=0;
y=0;
for i=1:n
    t(i)=deltaTdata(i,1);
    x(i,1)=1;
    x(i,2)=t(i)-1900;
    y(i,1)=deltaTdata(i,2);
end

disp('--------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,k]=size(x);
DegreesOfFreedom=n-k;
varName(1,:)='(Intercept)';
varName(2,:)='year - 1900';
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

for i=1:n
    bestFitTunEq(i)=t(i);
    bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1);
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
    title(['Antarctic temperature anomaly, ACF of model residuals']);
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

for i=1:n
    bestFitTacf(i)=t(i);
    bestFitYacf(i)=betaGLSacf(1,1)+x(i,2)*betaGLSacf(2,1);
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

for i=1:n
    bestFitTdw(i)=t(i);
    bestFitYdw(i)=betaGLSdw(1,1)+x(i,2)*betaGLSdw(2,1);
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

for i=1:n
    bestFitTmdw(i)=t(i);
    bestFitYmdw(i)=betaGLSmdw(1,1)+x(i,2)*betaGLSmdw(2,1);
end


disp('========================================================');


% Plot regression lines:
plot(t,y,'k.-',bestFitTunEq,bestFitYunEq,'r*:',bestFitTacf,bestFitYacf,'gx-',bestFitTdw,bestFitYdw,'bo-',bestFitTmdw,bestFitYmdw,'m+-');
legend('observed','OLS fit',['GLS ACF lag-1 AR1(',num2str(R(2),3),')'],['GLS Durbin-Watson AR1(',num2str(rhoDW,3),')'],['GLS Modified D-W AR1(',num2str(rhoDWM,3),')'],'Location','northwest');
title([num2str(n),' measured data points, and best fit lines, Antarctic temperature anomaly']);
xlabel('Year');
ylabel('Temperature anomaly, K');
axis([1850 2000 -1.5 1.5]);


