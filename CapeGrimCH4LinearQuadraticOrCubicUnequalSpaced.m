%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads Cape Grim CH4 data (unequally spaced at first, then does OLS and GLS for %
% linear (y=b0+b1.t), quadratic (y=b0+b1.t+b2.t^2) or cubic (y=b0+b1.t+b2.t^2+b3.t^3) models. %
% This code is for Matlab or Octave (no special libraries required).                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set one of these true, or otherwise both false for a linear model
doQuatratic=false;
doCubic=false;

% This file contains the monthly mean baseline (background) methane (CH4) concentration
% measured at the Cape Grim Baseline Air Pollution Station (CGBAPS), 
% located near the north-west tip of Tasmania, Australia.
% CGBAPS is funded and managed by the Australian Bureau of Meteorology, and
% the scientific program is jointly supervised with CSIRO Oceans & Atmosphere.

% The three columns of this csv file are: Date (decimal years), CH4 (ppb), Std Dev (ppb)
methane=csvread('CapeGrimMethane.csv');

[n,~]=size(methane);

%%%%%%%%%%%%%%%%%%%%%
% Assemble the data %
%%%%%%%%%%%%%%%%%%%%%
m=0;
t=0;  % sample coordinates
x=0;  % OLS/GLS matrix
y=0;  % correlated data
ds=0; % dataset
for i=1:n
    dateDecimal=methane(i,1);
    CH4ppb=methane(i,2);
        t(i)=dateDecimal;
        x(i,1)=1;
        x(i,2)=t(i)-2000;
        ds(i)=1; % all the one dataset
        if doQuatratic
            x(i,3)=x(i,2)*x(i,2);
        elseif doCubic
            x(i,3)=x(i,2)*x(i,2);
            x(i,4)=x(i,2)*x(i,2)*x(i,2);
        end
        y(i,1)=CH4ppb;
end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n,k]=size(x);
    DegreesOfFreedom=n-k;
    varName(1,:)='(Intercept)';
    varName(2,:)='year-2000  ';
    if doQuatratic
        varName(3,:)='(yr-2000)^2';
    elseif doCubic
        varName(3,:)='(yr-2000)^2';
        varName(4,:)='(yr-2000)^3';
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
    %disp('t-value = beta / std.error');
    disp('                 Estimate     Std. Error    t value');
    for m=1:k
        fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaOLS(m,1),stdErrorOLS(m,1),tValue(m,1));
    end
    disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k),' degrees of freedom']);

    
    bestFitYunEq=0;
    upperCIunEq=0;
    lowerCIunEq=0;
    CIunEq=0;
    if doCubic
        for i=1:n
            bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1)+x(i,3)*betaOLS(3,1)+x(i,4)*betaOLS(4,1);
            upperCIunEq(i)=bestFitYunEq(i)+1.96*residStdErrOLS;
            lowerCIunEq(i)=bestFitYunEq(i)-1.96*residStdErrOLS;
            CIunEq(i)=1.96*residStdErrOLS;
        end
    elseif doQuatratic
        for i=1:n
            bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1)+x(i,3)*betaOLS(3,1);
            upperCIunEq(i)=bestFitYunEq(i)+1.96*residStdErrOLS;
            lowerCIunEq(i)=bestFitYunEq(i)-1.96*residStdErrOLS;
            CIunEq(i)=1.96*residStdErrOLS;
        end
    else
        for i=1:n
            bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1);
            upperCIunEq(i)=bestFitYunEq(i)+1.96*residStdErrOLS;
            lowerCIunEq(i)=bestFitYunEq(i)-1.96*residStdErrOLS;
            CIunEq(i)=1.96*residStdErrOLS;
        end
    end

    
    disp('--------------------------------------------------------');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First find kd = rd/ro where rd = effective distance %
    % between the datasets, equation (4.7).               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if m>0
        kd=-log(1-sumSameTrds/(2*m*variance));
    else
        % no zero distance pairs - set kd =0
        kd=0;
    end
    disp(['kd = rd/ro = ',num2str(kd),' = (dataset distance)/ro - see (4.7)']);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find the sum of squares of residual differences from the closest %
    % residual (by time difference) in the SAME dataset, UNWEIGHTED        %
    % and divided by 2, to give RHS of (4.8).                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gammaRHSidw=0;
    rClosest=0;
    for i=1:n
        rMin=999999999;
        for j=1:n
            r=abs(x(i,2)-x(j,2));
            if r>0 && r<rMin && ds(i)==ds(j)
                % closest so far from this dataset
                rMin=r;
                rClosest(i)=r;
                %thisGamma=((residuals(i,1)-residuals(j,1))^2)/r;
                thisGamma=((residuals(i,1)-residuals(j,1))^2);
            end
        end
        gammaRHSidw=gammaRHSidw+thisGamma;
        %disp(['rMin: ',num2str(rMin)]);
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
        for i=1:n
            %gammaLHSidw=gammaLHSidw+(1-exp(-sqrt((rClosest(i)/ro)^2+kd*kd)))/rClosest(i);
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt((rClosest(i)/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*variance;
    end
    rd=kd*ro;
    disp('Durbin-Watson equivalent estimate (4.8):');
    if kd==0
        disp(['RBF time ro = ',num2str(ro),',  rd = ',num2str(rd),', so equivalent to (4.5)']);
    else
        disp(['RBF time ro = ',num2str(ro),',  dataset rd = ',num2str(rd)]);
    end


    disp('--------------------------------------------------------');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),',  dataset rd = ',num2str(rd)]);
    disp(['Assembling S: ',num2str(n),' x ',num2str(n),' with ',num2str(m),' pairs of co-located points']);
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
    %disp(['Inverting S']);
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
        for i=1:n
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt((rClosest(i)/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*LHSmultiplier;
    end
    rd=kd*ro;
    disp('Modified Durbin-Watson equivalent estimate (4.9):');
    if kd==0
        disp(['RBF time ro = ',num2str(ro),',  rd = ',num2str(rd),', so equivalent to (4.6)']);
    else
        disp(['RBF time ro = ',num2str(ro),',  dataset rd = ',num2str(rd)]);
    end

    disp('--------------------------------------------------------');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),',  dataset rd = ',num2str(rd)]);
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
    %disp(['Inverting S']);
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
    %disp('t-value = beta / std.error');
    disp('                 Estimate     Std. Error    t value');
    for m=1:k
        fprintf('%s  %12.6f  %12.6f  %9.3f\n',varName(m,:),betaGLSmdw3(m,1),stdErrorMDW3(m,1),tValue(m,1));
    end
    disp(['Residual standard error: ',num2str(residStdErrMDW3),' on ',num2str(n-k),' degrees of freedom']);
    
    bestFitYmdw3=0;
    upperCImdw3=0;
    lowerCImdw3=0;
    CImdw3=0;
    if doCubic
        for i=1:n
            bestFitYmdw3(i)=betaGLSmdw3(1,1)+x(i,2)*betaGLSmdw3(2,1)+x(i,3)*betaGLSmdw3(3,1)+x(i,4)*betaGLSmdw3(4,1);
            upperCImdw3(i)=bestFitYmdw3(i)+1.96*residStdErrMDW3;
            lowerCImdw3(i)=bestFitYmdw3(i)-1.96*residStdErrMDW3;
            CImdw3(i)=1.96*residStdErrMDW3;
        end
    elseif doQuatratic
        for i=1:n
            bestFitYmdw3(i)=betaGLSmdw3(1,1)+x(i,2)*betaGLSmdw3(2,1)+x(i,3)*betaGLSmdw3(3,1);
            upperCImdw3(i)=bestFitYmdw3(i)+1.96*residStdErrMDW3;
            lowerCImdw3(i)=bestFitYmdw3(i)-1.96*residStdErrMDW3;
            CImdw3(i)=1.96*residStdErrMDW3;
        end
    else
        for i=1:n
            bestFitYmdw3(i)=betaGLSmdw3(1,1)+x(i,2)*betaGLSmdw3(2,1);
            upperCImdw3(i)=bestFitYmdw3(i)+1.96*residStdErrMDW3;
            lowerCImdw3(i)=bestFitYmdw3(i)-1.96*residStdErrMDW3;
            CImdw3(i)=1.96*residStdErrMDW3;
        end
    end
    
    
    disp('========================================================');
    
    % Plot regression lines:
    plot(t,y,'r-',t,upperCIunEq,'g:',t,bestFitYunEq,'g-',t,lowerCIunEq,'g:',t,upperCImdw3,'b:',t,bestFitYmdw3,'b-',t,lowerCImdw3,'b:');
    legend('observed','95% OLS CI','OLS fit','95% OLS CI','95% GLS CI','GLS RBF fit','95% GLS CI','location','northwest');
    xlabel('Date');
    ylabel('Background methane, ppb');
    if doCubic
        title([num2str(n),' data points, and best fit lines - cubic model']);
    elseif doQuatratic
        title([num2str(n),' data points, and best fit lines - quadratic model']);
    else
        title([num2str(n),' data points, and best fit lines - linear model']);
    end
