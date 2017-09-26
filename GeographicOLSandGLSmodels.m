%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads multipath fading data, then does OLS and GLS linear regression.    %
% GLS assumes residual correlation is by geographic exponential radial basis functions. % 
% This code is for Matlab or Octave (clear or set the usingOctave flag).                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usingOctave=true; % forces flushing of the display buffer in Octave to see progress (not for Matlab)

% Read data:
multipath=csvread('Multipath.csv');
% This file contains the station numbers with their latitude & longitude,
% observed 0.01% of worst month of average year radio signal fading distribution,
% and the 9 model variables used in Salamon, Hansen & Abbott, ISAP 2016
% "Radio link clear-air fading prediction from surface weather station data."
% Columns of the csv file are:
% STATION, lat, long, A0d01, logD, dN1, log1plusEp, dN010ERAI, v2, v1, NsA0d1pc, HL, logFplus6
%    1      2    3      4     5     6        7          8      9   10     11     12     13

earthRad=6371; % used in distance calcs

%var=5;  % log(D)
var=10; % v1

    [n,~]=size(multipath);

    %%%%%%%%%%%%%%%%%%%%%
    % Assemble the data %
    %%%%%%%%%%%%%%%%%%%%%
    lat=0;  % data point latitude
    long=0; % data point longitude
    t=0;  % parameter v1
    x=0;  % OLS/GLS matrix
    y=0;  % correlated data
    ds=0; % dataset
    for i=1:n
        lat(i)=multipath(i,2);
        long(i)=multipath(i,3);
        t(i)=multipath(i,var); % current parameter
        x(i,1)=1;
        x(i,2)=t(i);
        % add x(i,3) etc if required for a multiple parameter model
        ds(i)=1; % all the one dataset
        y(i,1)=multipath(i,4); % a0.01%
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n,k]=size(x);
    DegreesOfFreedom=n-k;
    varName(1,:)='(Intercept)';
    if var==5
        varName(2,:)='log D      ';
    elseif var==6
        varName(2,:)='dN1 92-93  ';
    elseif var==7
        varName(2,:)='log(1+|Ep|)';
    elseif var==8
        varName(2,:)='dN010ERAI  ';
    elseif var==9
        varName(2,:)='v2         ';
    elseif var==10
        varName(2,:)='v1         ';
    elseif var==11
        varName(2,:)='NsA0.1%    ';
    elseif var==12
        varName(2,:)='HL         ';
    else % if var==13
        varName(2,:)='log(f+6)   ';
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

    
    bestFitYunEq=0;
    upperCIunEq=0;
    lowerCIunEq=0;
    CIunEq=0;
    for i=1:n
        bestFitYunEq(i)=betaOLS(1,1)+x(i,2)*betaOLS(2,1);
        upperCIunEq(i)=bestFitYunEq(i)+1.96*residStdErrOLS;
        lowerCIunEq(i)=bestFitYunEq(i)-1.96*residStdErrOLS;
        CIunEq(i)=1.96*residStdErrOLS;
    end

    
    disp('--------------------------------------------------------');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First find kd = rd/ro where rd = effective distance %
    % between the datasets, equation (3.7).               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=0;
    sumSameTrds=0;
    for i=1:n
        for j=1:n
            if abs(lat(i)-lat(j))+abs(long(i)-long(j))<0.001 && i~=j
                % same location, different data point
                m=m+1;
                sumSameTrds=sumSameTrds+(residuals(i,1)-residuals(j,1))^2;
            end
        end
    end
    variance=residStdErrOLS^2;
    if m>0
        kd=-log(1-sumSameTrds/(2*m*variance));
    else
        % no zero distance pairs - set kd =0
        kd=0;
    end
    disp(['kd = rd/ro = ',num2str(kd),' = (data same locn dist)/ro - see (4.7)']);
    
    if usingOctave
        fflush(stdout);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find the sum of squares of residual differences from the closest     %
    % residual (great circle distance) and divided by 2, to give RHS of (3.8). %
    % If there are more than one at the same nearest LOCATION, take the mean.  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gammaRHSidw=0;
    rClosest=0;
    nSamePlace=0;
    nDiffPlace=0;
    for i=1:n
        rMin=999999999;
        for j=1:n
            % calculate great-circle distance by spherical geometry
            phiS=(pi/180)*lat(i);
            phiT=(pi/180)*lat(j);
            lambdaD=(pi/180)*(long(i)-long(j));
            X=real(acos(sin(phiT)*sin(phiS)+cos(phiT)*cos(phiS)*cos(lambdaD)));
            r=X*earthRad;
            if r>0 && r<rMin && ds(i)==ds(j)
                % closest so far
                rMin=r;
                rClosest(i)=r;
                jClosest(i)=j;
            end
        end
        % check if there are any other residuals at the same coordinates as jClosest(i)
        m=0;
        nearestResid=0;
        for j=1:n
            if abs(lat(j)-lat(jClosest(i)))+abs(long(j)-long(jClosest(i)))<0.001
                % same closest location -> include in mean of closest residuals
                m=m+1;
                nearestResid=nearestResid+residuals(j,1);
            end
        end
        nearestResid=nearestResid/m;
        thisGamma=((residuals(i,1)-nearestResid)^2);
        if m>1
            nSamePlace=nSamePlace+1;
        else
            nDiffPlace=nDiffPlace+1;
        end

        gammaRHSidw=gammaRHSidw+thisGamma;
    end
    gammaRHSidw=gammaRHSidw/2;
    disp([num2str(100*nDiffPlace/(nSamePlace+nDiffPlace)),'% of data points are at unique locations']);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find ro by incrementing it until LHS of (3.8) <= RHS of (3.8) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=0;
    gammaLHSidw=2*gammaRHSidw;
    while gammaLHSidw>gammaRHSidw
        ro=ro+0.1;
        gammaLHSidw=0;
        for i=1:n
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt((rClosest(i)/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*variance;
    end
    rd=kd*ro;
    disp('Durbin-Watson equivalent estimate (4.8):');
    disp(['RBF distance ro = ',num2str(ro),' km,  data same locn rd = ',num2str(rd),' km']);


    disp('--------------------------------------------------------');
    
    if usingOctave
        fflush(stdout);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),' km,  data same locn rd = ',num2str(rd),' km']);
    disp(['Assembling S: ',num2str(n),' x ',num2str(n),', ',num2str(nDiffPlace),' data points are unique locations']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            % calculate great-circle distance by spherical geometry
            phiS=(pi/180)*lat(i);
            phiT=(pi/180)*lat(j);
            lambdaD=(pi/180)*(long(i)-long(j));
            X=real(acos(sin(phiT)*sin(phiS)+cos(phiT)*cos(phiS)*cos(lambdaD)));
            r=X*earthRad;
            phi=exp(-sqrt((r/ro)^2+kd*kd));
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    %disp(['Inverting S']);
    if usingOctave
        fflush(stdout);
    end
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


    disp('--------------------------------------------------------');
    
    if usingOctave
        fflush(stdout);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now find ro by incrementing it until LHS of (3.9) <= RHS of (3.9) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=0;
    LHSmultiplier=variance*(2-exp(-SSPD/SSR));
    gammaLHSidw=2*gammaRHSidw;
    while gammaLHSidw>gammaRHSidw
        ro=ro+0.1;
        gammaLHSidw=0;
        for i=1:n
            gammaLHSidw=gammaLHSidw+(1-exp(-sqrt((rClosest(i)/ro)^2+kd*kd)));
        end
        gammaLHSidw=gammaLHSidw*LHSmultiplier;
    end
    rd=kd*ro;
    disp('Modified Durbin-Watson equivalent estimate (4.9):');
    disp(['RBF distance ro = ',num2str(ro),' km,  data same locn rd = ',num2str(rd),' km']);

    disp('--------------------------------------------------------');
    
    if usingOctave
        fflush(stdout);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['GLS using ro = ',num2str(ro),' km,  data same locn rd = ',num2str(rd),' km']);
    disp(['Assembling S: ',num2str(n),' x ',num2str(n),', ',num2str(nDiffPlace),' data points are unique locations']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            % calculate great-circle distance by spherical geometry
            phiS=(pi/180)*lat(i);
            phiT=(pi/180)*lat(j);
            lambdaD=(pi/180)*(long(i)-long(j));
            X=real(acos(sin(phiT)*sin(phiS)+cos(phiT)*cos(phiS)*cos(lambdaD)));
            r=X*earthRad;
            phi=exp(-sqrt((r/ro)^2+kd*kd));
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    %disp(['Inverting S']);
    if usingOctave
        fflush(stdout);
    end
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
    
    bestFitYmdw3=0;
    upperCImdw3=0;
    lowerCImdw3=0;
    CImdw3=0;
    for i=1:n
        bestFitYmdw3(i)=betaGLSmdw3(1,1)+x(i,2)*betaGLSmdw3(2,1);
        upperCImdw3(i)=bestFitYmdw3(i)+1.96*residStdErrMDW3;
        lowerCImdw3(i)=bestFitYmdw3(i)-1.96*residStdErrMDW3;
        CImdw3(i)=1.96*residStdErrMDW3;
    end
        
    disp('========================================================');
    
    % Plot regression lines:
    plot(t,y,'r.',t,upperCIunEq,'g:',t,bestFitYunEq,'g-',t,lowerCIunEq,'g:',t,upperCImdw3,'b:',t,bestFitYmdw3,'b-',t,lowerCImdw3,'b:');
    legend('observed','95% OLS CI','OLS fit','95% OLS CI','95% GLS CI','GLS RBF fit','95% GLS CI','Location','northwest');
    title([num2str(n),' data points, and best fit lines']);
    xlabel(varName(2,:));
    ylabel('A0.01%, dB');


