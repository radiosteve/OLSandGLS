%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads monthly global temperature anomaly data,       %
% then does OLS and GLS for a linear (y = b0 + b1.t) model.         %
% This code is for Matlab or Octave (no special libraries required).%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annualMeans=false;    % take annual means, otherwise monthly

% choose one of these [first 2 columns of both are year, anomaly (K)]

% NOAAGlobalTemp .asc file of the same name but '   ' replaced with ',' to make it csv
TempAnomData=csvread('aravg.mon.land_ocean.90S.90N.v4.0.1.201704.csv');
dataName='NOAAGlobal';

% HADCrut4 .txt file of the same name but '   ' and '/' replaced with ',' to make it csv
%TempAnomData=csvread('HadCRUT.4.5.0.0.monthly_ns_avg.csv');
%dataName='HADCrut4';




    %%%%%%%%%%%%%%%%%%%%%
    % Assemble the data %
    %%%%%%%%%%%%%%%%%%%%%
    [n,~]=size(TempAnomData);
    t=0;
    x=0;
    y=0;
    m=0;
    for i=1:n
        tYear=TempAnomData(i,1)+(TempAnomData(i,2)-0.5)/12;
        % analyse the 40 years 1977-2016
        if tYear>1977 && tYear<2017
            if annualMeans
                if TempAnomData(i,2)==1
                    m=m+1;
                    t(m)=TempAnomData(i,1)+0.5;
                    x(m,1)=1;
                    x(m,2)=t(m)-2000;
                    y(m,1)=TempAnomData(i,3)/12;
                else
                    y(m,1)=y(m,1)+TempAnomData(i,3)/12;
                end
            else
                m=m+1;
                t(m)=tYear;
                x(m,1)=1;
                x(m,2)=t(m)-2000;
                y(m,1)=TempAnomData(i,3);
            end
        end
    end
    % adjust to zero mean over the period analysed
    yMean=mean(y);
    y=y-yMean;


    disp('========================================================');
    disp([dataName,' Annual Mean of Monthly Temperature Anomaly Data 1977-2016']);
    disp('========================================================');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n,k]=size(x);
    DegreesOfFreedom=n-k;
    varName(1,:)='(Intercept)';
    varName(2,:)='year - 2000';
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

    bestFitTunEq=0;
    bestFitYunEq=0;
    upperCIunEq=0;
    lowerCIunEq=0;
    for i=1:n
        bestFitTunEq(i)=t(i);
        bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1);
        upperCIunEq(i)=bestFitYunEq(i)+1.96*residStdErrOLS;
        lowerCIunEq(i)=bestFitYunEq(i)-1.96*residStdErrOLS;
    end
    
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Durbin-Watson Statistic, rho estimate and improved estimate %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SSRFD=0;
    for i=2:n
        SSRFD=SSRFD+(residuals(i)-residuals(i-1))^2;
    end
    dw=SSRFD/SSR;

    % Durbin-Watson estimate
    rhoDW=1-0.5*dw;

    % Modified Durbin-Watson estimate
    rhoDWM=1-0.5*SSRFD/(SSR*(2-exp(-SSPD/SSR)));
    
    variance=var(residuals);
    disp(['Residual variance = ',num2str(variance),' deg^2']);

    disp('--------------------------------------------------------');


    rho=rhoDW;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using Durbin-Watson for AR1 rho = ',num2str(rho)]);
    if annualMeans
        ro=-1/log(rho);
    else
        ro=-1/(12*log(rho));
    end
    disp(['        or as radial basis function, ro = ',num2str(ro),' years']);
    Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
    Sinv(1,1)=1/(1-rho*rho);
    Sinv(n,n)=1/(1-rho*rho);
    for i=1:n-1
        Sinv(i,i+1)=-rho/(1-rho*rho);
        Sinv(i+1,i)=-rho/(1-rho*rho);
    end
    betaGLSdw=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLSdw;
    residStdErr=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
    stdErrorDW=0;
    tValue=0;
    normCov=inv(transpose(x)*Sinv*x);
    for m=1:k
        stdErrorDW(m,1)=residStdErr*sqrt(normCov(m,m));
        tValue(m,1)=betaGLSdw(m,1)/stdErrorDW(m,1);
    end
    disp('                 Estimate     Std. Error    t value');
    for m=1:k
        fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSdw(m,1),stdErrorDW(m,1),tValue(m,1));
    end
    disp(['Residual standard error: ',num2str(residStdErr),' on ',num2str(n-k),' degrees of freedom']);
    
    disp('--------------------------------------------------------');



    rho=rhoDWM;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using modified Durbin-Watson for AR1 rho = ',num2str(rho)]);
    if annualMeans
        ro=-1/log(rho);
    else
        ro=-1/(12*log(rho));
    end
    disp(['        or as radial basis function, ro = ',num2str(ro),' years']);
    Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
    Sinv(1,1)=1/(1-rho*rho);
    Sinv(n,n)=1/(1-rho*rho);
    for i=1:n-1
        Sinv(i,i+1)=-rho/(1-rho*rho);
        Sinv(i+1,i)=-rho/(1-rho*rho);
    end
    betaGLSmdw3=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLSmdw3;
    residStdErrMDW3=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
    stdErrorMDW3=0;
    tValue=0;
    normCov=inv(transpose(x)*Sinv*x);
    for m=1:k
        stdErrorMDW3(m,1)=residStdErrMDW3*sqrt(normCov(m,m));
        tValue(m,1)=betaGLSmdw3(m,1)/stdErrorMDW3(m,1);
    end
    %disp('t-value = beta / std.error');
    disp('                 Estimate     Std. Error    t value');
    for m=1:k
        fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSmdw3(m,1),stdErrorMDW3(m,1),tValue(m,1));
    end
    disp(['Residual standard error: ',num2str(residStdErrMDW3),' on ',num2str(n-k),' degrees of freedom']);
    
    bestFitTmdw3=0;
    bestFitYmdw3=0;
    upperCImdw3=0;
    lowerCImdw3=0;
    for i=1:n
        bestFitTmdw3(i)=t(i);
        bestFitYmdw3(i)=betaGLSmdw3(1,1)+x(i,2)*betaGLSmdw3(2,1);
        upperCImdw3(i)=bestFitYmdw3(i)+1.96*residStdErrMDW3;
        lowerCImdw3(i)=bestFitYmdw3(i)-1.96*residStdErrMDW3;
    end
    
    
    disp('========================================================');

    
    % Plot regression lines:
    plot(t,y,'k-',bestFitTunEq,upperCIunEq,'b:',bestFitTunEq,bestFitYunEq,'b--',bestFitTunEq,lowerCIunEq,'b:',bestFitTmdw3,upperCImdw3,'r--',bestFitTmdw3,bestFitYmdw3,'r-',bestFitTmdw3,lowerCImdw3,'r--');
    legend('observed','OLS upper 95%CI','OLS fit','OLS lower 95%CI','OLS upper 95%CI',['GLS AR1(',num2str(rhoDWM,3),')'],'OLS lower 95%CI','Location','northwest');
    title([num2str(n),' measured data points, and best fit lines, ',dataName,'Temp Anomaly, (for zero mean 1977-2016)']);
    xlabel('Year');
    ylabel('Global Temperature Anomaly, K');

