%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads monthly global temperature anomaly data, then does OLS and %
% GLS for both linear (y = b0 + b1.t) and quadratic (y = b0 + b1.t + b2.t^2).   %
% This code is for Matlab or Octave (no special libraries required).            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***** This version takes the mean of each 12 months of the year to get annual data *****

dataName='NOAAGlobalAndHADCrut4';




    %%%%%%%%%%%%%%%%%%%%%
    % Assemble the data %
    %%%%%%%%%%%%%%%%%%%%%
    % NOAAGlobalTemp .asc file of the same name but '   ' replaced with ',' to make it csv
    TempAnomData1=csvread('aravg.mon.land_ocean.90S.90N.v4.0.1.201704.csv');
    [n1,~]=size(TempAnomData1);
    TempAnomData1a=0;
    m=0;
    for i=1:n1
        if TempAnomData1(i,1)>1976 && TempAnomData1(i,1)<2017
            m=m+1;
            TempAnomData1a(m,1)=TempAnomData1(i,1);
            TempAnomData1a(m,2)=TempAnomData1(i,2);
            TempAnomData1a(m,3)=TempAnomData1(i,3);
        end
    end
    Temp1Mean=mean(TempAnomData1a);
    % Temp1Mean(1,3) is the mean NOAA anomaly over 1977 to 2016

    % HADCrut4 .txt file of the same name but '   ' and '/' replaced with ',' to make it csv
    TempAnomData2=csvread('HadCRUT.4.5.0.0.monthly_ns_avg.csv');
    [n2,~]=size(TempAnomData2);
    TempAnomData2a=0;
    m=0;
    for i=1:n2
        if TempAnomData2(i,1)>1976 && TempAnomData2(i,1)<2017
            m=m+1;
            TempAnomData2a(m,1)=TempAnomData2(i,1);
            TempAnomData2a(m,2)=TempAnomData2(i,2);
            TempAnomData2a(m,3)=TempAnomData2(i,3);
        end
    end
    Temp2Mean=mean(TempAnomData2a);
    % Temp2Mean(1,3) is the mean HADCrut4 anomaly over 1977 to 2016
    
    % Combine the datasets
    % columns will be year, month, dataset, temp anomaly
    TempAnomData=0;
    m=0;
    [n1,~]=size(TempAnomData1a);
    for i=1:n1
        m=m+1;
        TempAnomData(m,1)=TempAnomData1a(i,1);
        TempAnomData(m,2)=TempAnomData1a(i,2);
        TempAnomData(m,3)=1;
        % adjust for zero mean 1977-2016
        TempAnomData(m,4)=TempAnomData1a(i,3)-Temp1Mean(1,3);
    end
    [n2,~]=size(TempAnomData2a);
    for i=1:n2
        m=m+1;
        TempAnomData(m,1)=TempAnomData2a(i,1);
        TempAnomData(m,2)=TempAnomData2a(i,2);
        TempAnomData(m,3)=2;
        % adjust for zero mean 1977-2016
        TempAnomData(m,4)=TempAnomData2a(i,3)-Temp2Mean(1,3);
    end
    % sort in year, month, dataset order
    TempAnomData=sortrows(TempAnomData);


    [n,~]=size(TempAnomData);
    t=0;  % time in decimal years AD
    x=0;  % OLS/GLS matrix
    y=0;  % temp anomaly
    ds=0; % dataset
    m=0;
    for i=1:n
        tYear=TempAnomData(i,1)+(TempAnomData(i,2)-0.5)/12;
        m=m+1;
        t(m)=tYear;
        x(m,1)=1;
        x(m,2)=t(m)-2000;
        ds(m)=TempAnomData(i,3);
        y(m,1)=TempAnomData(i,4);
    end



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
    
    disp('--------------------------------------------------------');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First find kd = rd/ro where rd = effective distance in years %
    % between the datasets, equation (4.7).                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=0;
    sumSameTrds=0;
    for i=2:n
        if abs(t(i)-t(i-1))<0.01 && ds(i)~=ds(i-1)
            % same time, different dataset
            m=m+1;
            sumSameTrds=sumSameTrds+(residuals(i,1)-residuals(i-1,1))^2;
        end
    end
    variance=residStdErrOLS^2;
    kd=-log(1-sumSameTrds/(2*m*variance));
    disp(['kd = rd/ro = ',num2str(kd),' = (dataset distance)/ro - see (4.7)']);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find the sum of squares of residual differences from the closest %
    % residual (by time difference) in the SAME dataset, inverse distance  %
    % weighted and divided by 2, to give RHS of (4.8).                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gammaRHSidw=0;
    for i=3:n
        if ds(i)~=ds(i-2)
            % checking - shouldn't happen
            disp(['Dataset alternation fail at row ',num2str(i),' (missing data?)']);
        end
        % unweighted:
        gammaRHSidw=gammaRHSidw+(residuals(i,1)-residuals(i-2,1))^2;
    end
    gammaRHSidw=gammaRHSidw/2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find ro by incrementing it until LHS of (4.8) <= RHS of (4.8) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=0;
    gammaLHSidw=2*gammaRHSidw;
    while gammaLHSidw>gammaRHSidw
        ro=ro+0.01;
        gammaLHSidw=0;
        for i=3:n
            if ds(i)~=ds(i-2)
                % checking - shouldn't happen
                disp(['Dataset alternation fail at row ',num2str(i),' (missing data?)']);
            end
            % unweighted:
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt(((x(i,2)-x(i-2,2))/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*variance;
    end
    rd=kd*ro;
    disp('Durbin-Watson equivalent estimate (4.8):');
    disp(['RBF time ro = ',num2str(ro),' years,  dataset rd = ',num2str(rd),' years']);


    disp('--------------------------------------------------------');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),' years,  dataset rd = ',num2str(rd),' years']);
    disp(['Assembling S: ',num2str(n),' x ',num2str(n)]);
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            if ds(i)==ds(j)
                phi=exp(-sqrt(((t(j)-t(i))/ro)^2+0));
            else
                phi=exp(-sqrt(((t(j)-t(i))/ro)^2+kd*kd));
            end
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    disp(['Inverting S']);
    Sinv=inv(S);
    betaGLSmdw3=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLSmdw3;
    residStdErr=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k));
    stdErrorMDW3=0;
    tValue=0;
    normCov=inv(transpose(x)*Sinv*x);
    for m=1:k
        stdErrorMDW3(m,1)=residStdErr*sqrt(normCov(m,m));
        tValue(m,1)=betaGLSmdw3(m,1)/stdErrorMDW3(m,1);
    end
    disp('                 Estimate     Std. Error    t value');
    for m=1:k
        fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSmdw3(m,1),stdErrorMDW3(m,1),tValue(m,1));
    end
    disp(['Residual standard error: ',num2str(residStdErr),' on ',num2str(n-k),' degrees of freedom']);

    disp('--------------------------------------------------------');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find ro by incrementing it until LHS of (4.9) <= RHS of (4.9) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=0;
    LHSmultiplier=variance*(2-exp(-SSPD/SSR));
    gammaLHSidw=2*gammaRHSidw;
    while gammaLHSidw>gammaRHSidw
        ro=ro+0.01;
        gammaLHSidw=0;
        for i=3:n
            if ds(i)~=ds(i-2)
                % checking - shouldn't happen
                disp(['Dataset alternation fail at row ',num2str(i),' (missing data?)']);
            end
            % unweighted:
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt(((x(i,2)-x(i-2,2))/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*LHSmultiplier;
    end
    rd=kd*ro;
    disp('Modified Durbin-Watson equivalent estimate (4.9):');
    disp(['RBF time ro = ',num2str(ro),' years,  dataset rd = ',num2str(rd),' years']);

    disp('--------------------------------------------------------');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),' years,  dataset rd = ',num2str(rd),' years']);
    disp(['Assembling S: ',num2str(n),' x ',num2str(n)]);
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            if ds(i)==ds(j)
                phi=exp(-sqrt(((t(j)-t(i))/ro)^2+0));
            else
                phi=exp(-sqrt(((t(j)-t(i))/ro)^2+kd*kd));
            end
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    disp(['Inverting S']);
    Sinv=inv(S);
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
    legend('observed','OLS upper 95%CI','OLS fit','OLS lower 95%CI','OLS upper 95%CI',['GLS AR1 fit'],'OLS lower 95%CI','Location','northwest');
    title([num2str(n),' measured data points, and best fit lines, ',dataName,'Temp Anomaly, (for zero mean 1977-2016)']);
    xlabel('Year');
    ylabel('Global Temperature Anomaly, K');


