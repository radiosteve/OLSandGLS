%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads monthly global temperature anomaly data, then does OLS and GLS for  %
% linear (modelOrder=1: y = b0 + b1.t), quadratic (modelOrder=2: y = b0 + b1.t + b2.t^2) %
% and higher order polynomial trends, with 95% prediction intervals.                     %
% This code is for Matlab or Octave (no special libraries required).                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usingOctave=false;

annualMeans=true;  % take annual means, otherwise monthly

modelOrder=2;      % regression model order (max 9 with this code, but 6 is max for inv(X'X) with Cape Grim data)

plotACF=false;     % plot autocorrelation function of OLS residuals?

plotBestFit=3;     % plot best fit with CI of mean response and PI: 0=no, 1=OLS, 2=ACF, 3=DW, 4=ML, 5=TADW, 6=DW-TADW extrapolated

predictYears=0;    % extend plotting (whole number of years)


% first 3 columns of both datasets are year, month, anomaly (K)

% NOAAGlobalTemp .asc file of the same name but '   ' replaced with ',' to make it csv
TempAnomData1=csvread('aravg.mon.land_ocean.90S.90N.v4.0.1.201704.csv');

% HADCrut4 .txt file of the same name but '   ' and '/' replaced with ',' to make it csv
TempAnomData2=csvread('HadCRUT.4.5.0.0.monthly_ns_avg.csv');

dataName='meanNOAAGlobal-HADCrut4';

[n1,~]=size(TempAnomData1);
[n2,~]=size(TempAnomData2);

%%%%%%%%%%%%%%%%%%%%%
% Assemble the data %
%%%%%%%%%%%%%%%%%%%%%
% synchronise the two datasets
m1=0;
m2=0;
tYear1=TempAnomData1(m1+1,1)+(TempAnomData1(m1+1,2)-0.5)/12;
tYear2=TempAnomData2(m2+1,1)+(TempAnomData2(m2+1,2)-0.5)/12;
while tYear1>tYear2
    m2=m2+1;
    tYear2=TempAnomData2(m2+1,1)+(TempAnomData2(m2+1,2)-0.5)/12;
end
while tYear2>tYear1
    m1=m1+1;
    tYear1=TempAnomData1(m1+1,1)+(TempAnomData1(m1+1,2)-0.5)/12;
end
if n1-m1>n2-m2
    n=n2-m2;
else
    n=n1-m1;
end
    t=0;
    x=0;
    y1=0;
    y2=0;
    m=0;
    for i=1:n
        tYear1=TempAnomData1(i+m1,1)+(TempAnomData1(i+m1,2)-0.5)/12;
        tYear2=TempAnomData2(i+m2,1)+(TempAnomData2(i+m2,2)-0.5)/12;
        if abs(tYear1-tYear2)>0.001
            disp(['Time mismatch at ',num2str(tYear1),', ',num2str(tYear2)]);
            if usingOctave
                fflush(stdout);
            end
        end
        tYear=(tYear1+tYear2)/2;
        % analyse the 120 years 1897-2016
        if tYear>1897 && tYear<2017
            if annualMeans
                if TempAnomData1(i+m1,2)==1
                    m=m+1;
                    t(m)=TempAnomData1(i+m1,1)+0.5;
                    x(m,1)=1;
                    if modelOrder>0
                        x(m,2)=t(m);
                    end
                    y1(m,1)=TempAnomData1(i+m1,3)/12;
                    y2(m,1)=TempAnomData2(i+m2,3)/12;
                else
                    y1(m,1)=y1(m,1)+TempAnomData1(i+m1,3)/12;
                    y2(m,1)=y2(m,1)+TempAnomData2(i+m2,3)/12;
                end
            else
                predictYears=0; % prediction not implemented for monthly
                m=m+1;
                t(m)=tYear;
                x(m,1)=1;
                if modelOrder>0
                    x(m,2)=t(m);
                end
                y1(m,1)=TempAnomData1(i+m1,3);
                y2(m,1)=TempAnomData2(i+m2,3);
            end
        end
    end
    % adjust each dataset for zero mean y
    meanY1=mean(y1);
    meanY2=mean(y2);
    y1=y1-meanY1;
    y2=y2-meanY2;
    y=(y1+y2)/2;
    % adjust to zero mean time over the period analysed
    [n,~]=size(x);
    means=mean(x);
    meanT=means(2);
    if modelOrder>0
        for i=1:n
            x(i,2)=x(i,2)-meanT;
            if modelOrder>1
                for j=2:modelOrder
                    x(i,j+1)=x(i,2)^j;
                end
            end
        end
    end


disp('=========================================================================');
if annualMeans
    disp([dataName,' Annual Mean of Monthly Temperature Anomaly Data ',num2str(floor(t(1))),'-',num2str(floor(t(n)))]);
    dataName=[dataName,'Annual'];
else
    disp([dataName,' Monthly Temperature Anomaly Data ',num2str(floor(t(1))),'-',num2str(floor(t(n)))]);
    dataName=[dataName,'Monthly'];
end
if modelOrder>1
    disp('The 2nd and higher order expanatory variables are decorrelated w.r.t. previous variables');
end
disp('=========================================================================');

% expanded t() and x() for prediction
tp=t;
xp=x;
if predictYears>0
    for i=n+1:n+predictYears
        tp(i)=tp(i-1)+1;
        xp(i,1)=1;
        if modelOrder>0
            xp(i,2)=tp(i)-meanT;
            for j=2:modelOrder
                xp(i,j+1)=xp(i,2)^j;
            end
        end
    end
end

% decorrelate the explanatory variables for accurate prediction intervals
if modelOrder>1
    for k=2:modelOrder
        % remove any correlation between x(:,k+1) and x(:,1:k) by subtracting
        % an OLS estimate of x(:,k+1) in terms of x(:,1:k)
        X=x(:,1:k);
        Xp=xp(:,1:k);
        Y=x(:,k+1);
        Yp=xp(:,k+1);
        betaDecorr=inv((transpose(X)*X))*transpose(X)*Y;
        x(:,k+1)=Y-X*betaDecorr;
        xp(:,k+1)=Yp-Xp*betaDecorr;
    end
end

[n,k1]=size(x); % k1 = k + 1 = no. of variables + intercept
k=k1-1;
 
% A for SSRFD = z'A.z (see Durbin-Watson II 1951)
A=2*eye(n);
A(1,1)=1;
A(n,n)=1;
for i=2:n
    A(i-1,i)=-1;
    A(i,i-1)=-1;
end

% calculate dx = Xi - Xmean, for plotting prediction interval
meanX=mean(x);
meanX(1,1)=0; % except dx(intercept) = 1
dx=meanX*0;
for i=1:n
    dx(i,:)=x(i,:)-meanX;
end
dxp=dx;
if predictYears>0
    for i=n+1:n+predictYears
        dxp(i,:)=xp(i,:)-meanX;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varName(1,:)='(Intercept)';
if modelOrder>0
    varName(2,:)='year - mean';
end
if modelOrder>1
    for j=2:modelOrder
        varName(j+1,:)=['(yr-mean)^',num2str(j)];
    end
end
disp('OLS calculation by matrix algebra:');
invXdX=inv(transpose(x)*x);
betaOLS=invXdX*transpose(x)*y;
residuals=y-x*betaOLS;

SSR=transpose(residuals)*residuals;
residStdErrOLS=sqrt(transpose(residuals)*residuals/(n-k1));
stdErrorOLS=0;
tValue=0;
for m=1:k1
    stdErrorOLS(m,1)=residStdErrOLS*sqrt(invXdX(m,m));
    tValue(m,1)=betaOLS(m,1)/stdErrorOLS(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaOLS(m,1),stdErrorOLS(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residuals)*residuals)/n);
disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and plot sample autocorrelation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample autocorrelation function (provided in MatLab but not Octave)
Ymean=mean(residuals);
R=0;
Rall=0;
Lag=0;
zeroLine=0;
CIu=0;
CIl=0;
denom=0;
for i=1:n
    denom=denom+(residuals(i)-Ymean)^2;
end
RkSS=0;
for ka=0:n-1 %20
    numer=0;
    for na=1:n-ka
        numer=numer+(residuals(na)-Ymean)*(residuals(na+ka)-Ymean);
    end
    Rall(ka+1,1)=numer/denom;
    if ka<21
        R(ka+1,1)=numer/denom;
        Lag(ka+1)=ka;
        zeroLine(ka+1)=0;
        CIu(ka+1)=1.96/sqrt(n);
        CIl(ka+1)=-1.96/sqrt(n);
    end
    if ka>0
        RkSS=RkSS+Rall(ka+1,1)^2;
    end
end
if plotACF
    plot(Lag,R,'.-b',Lag,zeroLine,'-g',Lag,CIu,'-r',Lag,CIl,'-r');
    axis([0 20 -0.4 1]);
    xlabel('Lag');
    ylabel('Sample Autocorrelation');
    grid('on');
    title(['Cape Grim CO_2 annual data, ACF of order ',num2str(modelOrder),' model residuals']);
    csvwrite([dataName,'Order',num2str(modelOrder),'polynomialOLStrendACF.csv'],R);
end

SSR=transpose(residuals)*residuals;
SSRFD=0;
for i=2:n
    SSRFD=SSRFD+(residuals(i)-residuals(i-1))^2;
end
d=SSRFD/SSR;

% work out confidence limits of d
p=trace(A)-trace(transpose(x)*A*x*invXdX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*A*A*x*invXdX)+trace((transpose(x)*A*x*invXdX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
% approximate distribution of d by a beta disbn with same mean & variance (Durbin-Watson II 1951, III 1971)
ab=Ed*(4-Ed)/Vd-1;
a=ab*Ed/4;
b=ab-a;
d025=4*betainv(0.025,a,b);
d975=4*betainv(0.975,a,b);

if d<d025 || d>d975
    disp(['OLS residual Durbin-Watson: ',num2str(d),' [correl. error or bad model: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
else
    disp(['OLS residual Durbin-Watson: ',num2str(d),' [good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);

    bestFitUnEq=xp*betaOLS;
    tValue1=tinv(0.975,n-k-1);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),n-k-1);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrOLS^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorOLS(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitUnEq(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitUnEq(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitUnEq(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==1
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitUnEq(i,1)+predInt(i,1) || y1(i,1)<bestFitUnEq(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitUnEq(i,1)+predInt(i,1) || y2(i,1)<bestFitUnEq(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-','b.-',tp,bestFitUnEq,'r-',tp,bestFitUnEq+predInt,'g-',tp,bestFitUnEq-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (OLS)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [OLS residuals, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [OLS residuals, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitUnEq(i,1);               % predicted value
            CSVoutput(i,5)=bestFitUnEq(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitUnEq(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendOLS.csv'],CSVoutput);
    end

end

rho=R(2,1); % for GLS

disp('-------------------------------------------------------------------------');

if usingOctave
    fflush(stdout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sb scaling (similar Bartlett 1935 method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaling=sqrt((1+rho)/(1-rho));
disp(['Regression coeff std errors scaled by sqrt[(1+r1)/(1-r1)] = ',num2str(scaling)]);
stdErrorScaled=0;
tValue=0;
for m=1:k1
    stdErrorScaled(m,1)=scaling*residStdErrOLS*sqrt(invXdX(m,m));
    tValue(m,1)=betaOLS(m,1)/stdErrorScaled(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaOLS(m,1),stdErrorScaled(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);
disp('No transformed model Durbin-Watson test available as this method uses OLS');

disp('-------------------------------------------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['GLS using ACF lag 1 for AR1 rho = ',num2str(rho),' (lag 2, 3: ',num2str(R(3,1)),', ',num2str(R(4,1)),')']);

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

invXdSiX=inv(transpose(x)*Sinv*x);
betaGLSacf=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSacf;
residStdErrACF=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k1));
stdErrorACF=0;
tValue=0;
for m=1:k1
    stdErrorACF(m,1)=residStdErrACF*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLSacf(m,1)/stdErrorACF(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaGLSacf(m,1),stdErrorACF(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrACF),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

bestFitACF=x*betaGLSacf;

P=Sinv^0.5;
residualsEq=P*residualsGLS;
SSReq=transpose(residualsEq)*residualsEq;
SSRFDeq=0;
for i=2:n
    SSRFDeq=SSRFDeq+(residualsEq(i)-residualsEq(i-1))^2;
end
d=SSRFDeq/SSReq;

% work out confidence limits of d
p=trace(A)-trace(transpose(x)*transpose(P)*A*P*x*invXdSiX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*transpose(P)*A*A*P*x*invXdSiX)+trace((transpose(x)*transpose(P)*A*P*x*invXdSiX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
% approximate distribution of d by a beta disbn with same mean & variance (Durbin-Watson II 1951, III 1971)
ab=Ed*(4-Ed)/Vd-1;
a=ab*Ed/4;
b=ab-a;
d025=4*betainv(0.025,a,b);
d975=4*betainv(0.975,a,b);

if d<d025 || d>d975
    disp(['Transformed residual Durbin-Watson: ',num2str(d),' [bad model or eq: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
else
    disp(['Transformed residual Durbin-Watson: ',num2str(d),'   [good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);

    bestFitACF=xp*betaGLSacf;
    tValue1=tinv(0.975,n-k-1);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),n-k-1);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrACF^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorACF(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitACF(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitACF(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitACF(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==2
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitACF(i,1)+predInt(i,1) || y1(i,1)<bestFitACF(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitACF(i,1)+predInt(i,1) || y2(i,1)<bestFitACF(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitACF,'r-',tp,bestFitACF+predInt,'g-',tp,bestFitACF-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (ACF)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitACF(i,1);               % predicted value
            CSVoutput(i,5)=bestFitACF(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitACF(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendACF.csv'],CSVoutput);
    end
        
end


disp('-------------------------------------------------------------------------');

if usingOctave
    fflush(stdout);
end


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

% Thiel-Nagar (1961) estimate of rho
%rhoTN=((n^2)*(1-0.5*dw)+k1^2)/(n^2-k1^2);

% simple unbiased
%rhoUB=1-0.5*dw*(1+1/n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=rhoDW;
disp(['GLS using Durbin-Watson for AR1 rho = 1-d/2 = ',num2str(rho)]);
%disp(['GLS using Durbin-Watson for AR1 rho = ',num2str(rho),' (c.f. Thiel-Nagar: ',num2str(rhoTN),')']);
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
invXdSiX=inv(transpose(x)*Sinv*x);
betaGLSdw=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSdw;
residStdErrDW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k1));
stdErrorDW=0;
tValue=0;
for m=1:k1
    stdErrorDW(m,1)=residStdErrDW*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLSdw(m,1)/stdErrorDW(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaGLSdw(m,1),stdErrorDW(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrDW),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

P=Sinv^0.5;
residualsEq=P*residualsGLS;
SSReq=transpose(residualsEq)*residualsEq;
SSRFDeq=0;
for i=2:n
    SSRFDeq=SSRFDeq+(residualsEq(i)-residualsEq(i-1))^2;
end
d=SSRFDeq/SSReq;

% work out confidence limits of d
p=trace(A)-trace(transpose(x)*transpose(P)*A*P*x*invXdSiX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*transpose(P)*A*A*P*x*invXdSiX)+trace((transpose(x)*transpose(P)*A*P*x*invXdSiX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
% approximate distribution of d by a beta disbn with same mean & variance (Durbin-Watson II 1951, III 1971)
ab=Ed*(4-Ed)/Vd-1;
a=ab*Ed/4;
b=ab-a;
d025=4*betainv(0.025,a,b);
d975=4*betainv(0.975,a,b);

if d<d025 || d>d975
    disp(['Transformed residual Durbin-Watson: ',num2str(d),' [bad model or eq: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
else
    disp(['Transformed residual Durbin-Watson: ',num2str(d),'   [good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);

    bestFitDW=xp*betaGLSdw;
    tValue1=tinv(0.975,n-k-1);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),n-k-1);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrDW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorDW(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitDW(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitDW(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitDW(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==3
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitDW(i,1)+predInt(i,1) || y1(i,1)<bestFitDW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitDW(i,1)+predInt(i,1) || y2(i,1)<bestFitDW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitDW,'r-',tp,bestFitDW+predInt,'g-',tp,bestFitDW-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (DW)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitDW(i,1);               % predicted value
            CSVoutput(i,5)=bestFitDW(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitDW(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendDW.csv'],CSVoutput);
    end
        
end

disp('-------------------------------------------------------------------------');

if usingOctave
    fflush(stdout);
end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % max likelihood estimate (Shumway 4th ed eq 3.116 p118) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x1=0;
        for na=1:n-1
            x1=x1+residuals(na);
        end
        x1=x1/(n-1);
        x2=0;
        for na=2:n
            x2=x2+residuals(na);
        end
        x2=x2/(n-1);
        
        denom=0;
        for i=2:n
            denom=denom+(residuals(i-1)-x1)^2;
        end
        numer=0;
        for i=2:n
            numer=numer+(residuals(i)-x2)*(residuals(i-1)-x1);
        end
        rhoi=numer/denom;
    
        notConverged=true;
    
        while notConverged
        
            Sinv=eye(n)*(1+rhoi*rhoi)/(1-rhoi*rhoi);
            Sinv(1,1)=1/(1-rhoi*rhoi);
            Sinv(n,n)=1/(1-rhoi*rhoi);
            for i=1:n-1
                Sinv(i,i+1)=-rhoi/(1-rhoi*rhoi);
                Sinv(i+1,i)=-rhoi/(1-rhoi*rhoi);
            end
            betaGLSi=inv(transpose(x)*Sinv*x)*transpose(x)*Sinv*y;
            residualsGLSi=y-x*betaGLSi;
        
            % max likelihood estimate (Shumway 4th ed eq 3.116 p118)
            
            x1=0;
            for na=1:n-1
                x1=x1+residualsGLSi(na);
            end
            x1=x1/(n-1);
            x2=0;
            for na=2:n
                x2=x2+residualsGLSi(na);
            end
            x2=x2/(n-1);
            
            denom=0;
            for i=2:n
                denom=denom+(residualsGLSi(i)-x1)^2;
            end
            numer=0;
            for na=2:n
                numer=numer+(residualsGLSi(na)-x2)*(residualsGLSi(na-1)-x1);
            end
            rhon=numer/denom;
    
            if abs(rhoi-rhon)<0.000001
                notConverged=false;
            end
        
            rhoi=rhon;
        
        end
        
        rhoMLE=rhoi;


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=rhoMLE;
disp(['GLS using maximum likelihood for AR1 rho = ',num2str(rho)]);
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
invXdSiX=inv(transpose(x)*Sinv*x);
betaGLSmle=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSmle;
residStdErrMLE=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k1));
stdErrorMLE=0;
tValue=0;
for m=1:k1
    stdErrorMLE(m,1)=residStdErrMLE*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLSmle(m,1)/stdErrorMLE(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaGLSmle(m,1),stdErrorMLE(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrMLE),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);


P=Sinv^0.5;
residualsEq=P*residualsGLS;
SSReq=transpose(residualsEq)*residualsEq;
SSRFDeq=0;
for i=2:n
    SSRFDeq=SSRFDeq+(residualsEq(i)-residualsEq(i-1))^2;
end
d=SSRFDeq/SSReq;

% work out confidence limits of d
p=trace(A)-trace(transpose(x)*transpose(P)*A*P*x*invXdSiX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*transpose(P)*A*A*P*x*invXdSiX)+trace((transpose(x)*transpose(P)*A*P*x*invXdSiX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
% approximate distribution of d by a beta disbn with same mean & variance (Durbin-Watson II 1951, III 1971)
ab=Ed*(4-Ed)/Vd-1;
a=ab*Ed/4;
b=ab-a;
d025=4*betainv(0.025,a,b);
d975=4*betainv(0.975,a,b);

if d<d025 || d>d975
    disp(['Transformed residual Durbin-Watson: ',num2str(d),' [bad model or eq: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
else
    disp(['Transformed residual Durbin-Watson: ',num2str(d),'   [good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);

    bestFitMLE=xp*betaGLSmle;
    tValue1=tinv(0.975,n-k-1);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),n-k-1);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrMLE^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorMLE(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitMLE(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitMLE(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitMLE(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==4
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitMLE(i,1)+predInt(i,1) || y1(i,1)<bestFitMLE(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitMLE(i,1)+predInt(i,1) || y2(i,1)<bestFitMLE(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitMLE,'r-',tp,bestFitMLE+predInt,'g-',tp,bestFitMLE-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (MLE)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitMLE(i,1);               % predicted value
            CSVoutput(i,5)=bestFitMLE(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitMLE(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendMLE.csv'],CSVoutput);
    end
    
end


disp('-------------------------------------------------------------------------');


if usingOctave
    fflush(stdout);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tanh adjusted Durbin-Watson %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=(transpose(residuals)*A*residuals)/(transpose(residuals)*residuals);
% mean Ed and variance Vd of d
invXdX=inv(transpose(x)*x);
p=trace(A)-trace(transpose(x)*A*x*invXdX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*A*A*x*invXdX)+trace((transpose(x)*A*x*invXdX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
atanhRho=(2/(n-k1-3))*sqrt((n-k1+3)/Vd)*(atanh(1-d/2)-atanh(1-Ed/2));
rho=real(tanh(atanhRho));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['tanh adjusted Durbin-Watson, rho = ',num2str(rho)]);
Sinv=eye(n)*(1+rho*rho)/(1-rho*rho);
Sinv(1,1)=1/(1-rho*rho);
Sinv(n,n)=1/(1-rho*rho);
for i=1:n-1
    Sinv(i,i+1)=-rho/(1-rho*rho);
    Sinv(i+1,i)=-rho/(1-rho*rho);
end
invXdSiX=inv(transpose(x)*Sinv*x);
betaGLStadw=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLStadw;
residStdErrTADW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(n-k1));
stdErrorTADW=0;
tValue=0;
for m=1:k1
    stdErrorTADW(m,1)=residStdErrTADW*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLStadw(m,1)/stdErrorTADW(m,1);
end
disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaGLStadw(m,1),stdErrorTADW(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrTADW),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

P=Sinv^0.5;
residualsEq=P*residualsGLS;
SSReq=transpose(residualsEq)*residualsEq;
SSRFDeq=0;
for i=2:n
    SSRFDeq=SSRFDeq+(residualsEq(i)-residualsEq(i-1))^2;
end
d=SSRFDeq/SSReq;

% work out confidence limits of d
p=trace(A)-trace(transpose(x)*transpose(P)*A*P*x*invXdSiX);
Ed=p/(n-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(x)*transpose(P)*A*A*P*x*invXdSiX)+trace((transpose(x)*transpose(P)*A*P*x*invXdSiX)^2);
Vd=2*(q-p*Ed)/((n-k-1)*(n-k+1)); % variance of d (exact)
% approximate distribution of d by a beta disbn with same mean & variance (Durbin-Watson II 1951, III 1971)
ab=Ed*(4-Ed)/Vd-1;
a=ab*Ed/4;
b=ab-a;
d025=4*betainv(0.025,a,b);
d975=4*betainv(0.975,a,b);

if d<d025 || d>d975
    disp(['Transformed residual Durbin-Watson: ',num2str(d),' [bad model or eq: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
else
    disp(['Transformed residual Durbin-Watson: ',num2str(d),'   [good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);

    bestFitTADW=xp*betaGLStadw;
    tValue1=tinv(0.975,n-k-1);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),n-k-1);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrTADW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorTADW(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitTADW(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitTADW(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitTADW(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==5
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitTADW(i,1)+predInt(i,1) || y1(i,1)<bestFitTADW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitTADW(i,1)+predInt(i,1) || y2(i,1)<bestFitTADW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitTADW,'r-',tp,bestFitTADW+predInt,'g-',tp,bestFitTADW-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (TADW)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitTADW(i,1);               % predicted value
            CSVoutput(i,5)=bestFitTADW(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitTADW(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendTADW.csv'],CSVoutput);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply t and Sb extrapolation since d test ok %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('----- Regression coeff std errors and t extrapolated: 2*TADW-DW -----');
    stdErrorExtrap=0;
    tValue=0;
    for m=1:k+1
        stdErrorExtrap(m,1)=2*stdErrorTADW(m,1)-stdErrorDW(m,1);
        tValue(m,1)=2*(betaGLStadw(m,1)/stdErrorTADW(m,1))-(betaGLSdw(m,1)/stdErrorDW(m,1));
    end
    disp('                 Estimate     Std. Error    t value    Pr(>|t|)');
    for m=1:k1
        fprintf('%s  %12.6f  %12.6f  %9.3f  %12.6f\n',varName(m,:),betaGLStadw(m,1),stdErrorExtrap(m,1),tValue(m,1),2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
    end

    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n+predictYears
        predInt(i,1)=residStdErrTADW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorExtrap(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(n+predictYears)),': ',num2str(bestFitTADW(n+predictYears,1)),' deg; 95% prediction interval: ',num2str(bestFitTADW(n+predictYears,1)-predInt(n+predictYears,1)),', ',num2str(bestFitTADW(n+predictYears,1)+predInt(n+predictYears,1))]);
    end

    if plotBestFit==6
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitTADW(i,1)+predInt(i,1) || y1(i,1)<bestFitTADW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitTADW(i,1)+predInt(i,1) || y2(i,1)<bestFitTADW(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitTADW,'r-',tp,bestFitTADW+predInt,'g-',tp,bestFitTADW-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (DW-TADW extrap)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:n+predictYears
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitTADW(i,1);               % predicted value
            CSVoutput(i,5)=bestFitTADW(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitTADW(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendDWTADWextrap.csv'],CSVoutput);
    end

end

disp('=========================================================================');

