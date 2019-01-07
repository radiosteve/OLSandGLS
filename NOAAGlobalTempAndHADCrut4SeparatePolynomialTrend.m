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

plotBestFit=2;     % plot best fit with CI of mean response and PI: 0=no, 1=OLS, 2=DW, 3=TADW, 4=DW-TADW extrapolated, 5=Semivariogram

predictYears=0;    % extend plotting (whole number of years)


% first 3 columns of both datasets are year, month, anomaly (K)

% NOAAGlobalTemp .asc file of the same name but '   ' replaced with ',' to make it csv
TempAnomData1=csvread('aravg.mon.land_ocean.90S.90N.v4.0.1.201704.csv');

% HADCrut4 .txt file of the same name but '   ' and '/' replaced with ',' to make it csv
TempAnomData2=csvread('HadCRUT.4.5.0.0.monthly_ns_avg.csv');

dataName='sepRdRssDfSvgNOAAGlobal-HADCrut4';

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
    t=0;  % time for plotting
    x=0;  % X for combined data
    y1=0; % observed dataset 1 for plotting
    y2=0; % observed dataset 2 for plotting
    y=0;  % combined observations
    m=0;
    ds=0; % dataset number
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
                    x(2*m-1,1)=1;
                    x(2*m,1)=1;
                    if modelOrder>0
                        x(2*m-1,2)=t(m);
                        x(2*m,2)=t(m);
                    end
                    y1(m,1)=TempAnomData1(i+m1,3)/12;
                    y2(m,1)=TempAnomData2(i+m2,3)/12;
                    y(2*m-1,1)=y1(m,1);
                    y(2*m,1)=y2(m,1);
                    ds(2*m-1)=1;
                    ds(2*m)=2;
                else
                    y1(m,1)=y1(m,1)+TempAnomData1(i+m1,3)/12;
                    y2(m,1)=y2(m,1)+TempAnomData2(i+m2,3)/12;
                    y(2*m-1,1)=y(2*m-1,1)+TempAnomData1(i+m1,3)/12;
                    y(2*m,1)=y(2*m,1)+TempAnomData2(i+m2,3)/12;
                end
            else
                predictYears=0; % prediction not implemented for monthly
                m=m+1;
                t(m)=tYear;
                x(2*m-1,1)=1;
                x(2*m,1)=1;
                if modelOrder>0
                    x(2*m-1,2)=t(m);
                    x(2*m,2)=t(m);
                end
                y1(m,1)=TempAnomData1(i+m1,3);
                y2(m,1)=TempAnomData2(i+m2,3);
                y(2*m-1,1)=y1(m,1);
                y(2*m,1)=y2(m,1);
                ds(2*m-1)=1;
                ds(2*m)=2;
            end
        end
    end
    % adjust each dataset for zero mean y
    meanY1=mean(y1);
    meanY2=mean(y2);
    y1=y1-meanY1;
    y2=y2-meanY2;
    [n,~]=size(x);
    for i=1:n
        if ds(i)<1.5
            %disp(['y1: ',num2str(y(i,1)),' -> ',num2str(y(i,1)-meanY1)]);
            y(i,1)=y(i,1)-meanY1;
        else
            %disp(['y2: ',num2str(y(i,1)),' -> ',num2str(y(i,1)-meanY2)]);
            y(i,1)=y(i,1)-meanY2;
        end
    end
    
    % adjust to zero mean time over the period analysed
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
    disp([dataName,' Annual Mean of Monthly Temperature Anomaly Data ',num2str(floor(t(1))),'-',num2str(floor(t(n/2)))]);
    dataName=[dataName,'Annual'];
else
    disp([dataName,' Monthly Temperature Anomaly Data ',num2str(floor(t(1))),'-',num2str(floor(t(n/2)))]);
    dataName=[dataName,'Monthly'];
end
if modelOrder>1
    disp('The 2nd and higher order expanatory variables are decorrelated w.r.t. previous variables');
end
disp('=========================================================================');

% expanded t() and x() for prediction
tp=t;
[~,np]=size(t);
xp=0;
for i=1:np
    xp(i,1)=1;
    if modelOrder>0
        xp(i,2)=tp(i)-meanT;
    end
    if modelOrder>1
        for j=2:modelOrder
            xp(i,j+1)=xp(i,2)^j;
        end
    end
end
if predictYears>0
    for i=np+1:np+predictYears
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
[~,np]=size(tp);

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
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make list of "adjacent" points for "forward differences" for generalised Durbin-Watson %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making progressive list of nearest neighbours ...');
if usingOctave
    fflush(stdout);
end
% Find starting point with the greatest sum of distances to all other points
rSumMax=0;
nSumMax=0;
for i=1:n
    rSumThisOne=0;
    for j=1:n
        %r=sqrt((t(i)-t(j))^2+(u(i)-u(j))^2+(v(i)-v(j))^2); % 3D
        r=abs(x(i,2)-x(j,2)); % 1D
        rSumThisOne=rSumThisOne+r;
    end
    if rSumThisOne>rSumMax
        rSumMax=rSumThisOne;
        nSumMax=i;
    end
end
%disp(['Starting point is no. ',num2str(nSumMax),' at t = ',num2str(x(nSumMax,2))]);
%if usingOctave
%    fflush(stdout);
%end

% also generate reduction matrix to pre-multiply x, y or e to average co-located points
R=zeros(1,n); % first row
Rrow=1;
nThisRow=1;
R(Rrow,nSumMax)=1;

adjacentPoint=zeros(1,n);
adjacentDist=zeros(1,n);
%adjacentDist2=zeros(1,n);
thisPoint=nSumMax; % starting point
startPoint=nSumMax;
% find next point
rMin=999999999;
for j=1:n
    %r=sqrt((t(thisPoint)-t(j))^2+(u(thisPoint)-u(j))^2+(v(thisPoint)-v(j))^2); % 3D
    r=abs(x(thisPoint,2)-x(j,2)); % 1D
    %if r>0.001 && r<rMin
    if r<rMin && j~=thisPoint
        % closest so far
        rMin=r;
        adjacentPoint(thisPoint)=j;
        adjacentDist(thisPoint)=r;
        %adjacentDist2(thisPoint)=r*r;
    end
end
% and then continue until no more adjacent but different location points found
while adjacentPoint(thisPoint)>0
    thisPoint=adjacentPoint(thisPoint);
    if rMin>0
        if nThisRow>1
            % multiple entries on last row so divide it by nThisRow
            for i=1:n
                R(Rrow,i)=R(Rrow,i)/nThisRow;
            end
        end
        % start next row of R
        Rrow=Rrow+1;
        nThisRow=1;
    else
        % another entry at the same location; same row of R
        nThisRow=nThisRow+1;
    end
    R(Rrow,thisPoint)=1;
    adjacentPoint(thisPoint)=0;
    rMin=999999999;
    for j=1:n
        if adjacentPoint(j)==0
            %r=sqrt((t(thisPoint)-t(j))^2+(u(thisPoint)-u(j))^2+(v(thisPoint)-v(j))^2); % 3D
            r=abs(x(thisPoint,2)-x(j,2)); % 1D
            %if r>0.001 && r<rMin
            if r<rMin && j~=thisPoint
                % closest so far
                rMin=r;
                adjacentPoint(thisPoint)=j;
                adjacentDist(thisPoint)=r;
                %adjacentDist2(thisPoint)=r*r;
            end
        end
    end
end
% make sum of last row of R = 1
if nThisRow>1
    % multiple entries on last row so divide it by nThisRow
    for i=1:n
        R(Rrow,i)=R(Rrow,i)/nThisRow;
    end
end

%==== NEW VERSION ====
% use premultiplication matrix R to merge co-located points
% and sort in nearest new neighbour order
[nr,~]=size(R);
effDist=sum(adjacentDist)/(nr-1); % the last point on the path has adjacentDist = 0
disp(['Starting with no ',num2str(nSumMax),', mean non-zero distance = ',num2str(effDist),' years']);
disp('-------------------------------------------------------------------------');
if usingOctave
    fflush(stdout);
end
rx=R*x;
% Construct matrix A for Durbin-Watson d = (e'Ae)/(e'e)
A=2*eye(nr);
A(1,1)=1;
A(nr,nr)=1;
for i=1:nr-1
    A(i,i+1)=-1;
    A(i+1,i)=-1;
end


% calculate dx = Xi - Xmean, for plotting prediction interval
meanX=mean(x);
meanX(1,1)=0; % except dx(intercept) = 1
dx=meanX*0;
for i=1:n
    dx(i,:)=x(i,:)-meanX;
end
dxp=meanX*0;
for i=1:np
    dxp(i,:)=xp(i,:)-meanX;
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
residualsOLS=y-x*betaOLS;

SSR=transpose(residualsOLS)*residualsOLS;
residStdErrOLS=sqrt(transpose(residualsOLS)*residualsOLS/(n-k-1));
stdErrorOLS=0;
tValue=0;
for m=1:k1
    stdErrorOLS(m,1)=residStdErrOLS*sqrt(invXdX(m,m));
    tValue(m,1)=betaOLS(m,1)/stdErrorOLS(m,1);
end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k+1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaOLS(m,1),stdErrorOLS(m,1),tValue(m,1),n-k-1,2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsOLS)*residualsOLS)/n);
disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k-1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

SSRFD=(transpose(R*residualsOLS)*A*(R*residualsOLS));
SSR=transpose(R*residualsOLS)*(R*residualsOLS);
d=SSRFD/SSR;
% rho from Durbin-Watson
rhoDW=1-0.5*d;

% work out confidence limits of d
Rx=R*x;
invRXdRX=inv(transpose(Rx)*Rx);
p=trace(A)-trace(transpose(Rx)*A*Rx*invRXdRX);
Ed=p/(nr-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(Rx)*A*A*Rx*invRXdRX)+trace((transpose(Rx)*A*Rx*invRXdRX)^2);
Vd=2*(q-p*Ed)/((nr-k-1)*(nr-k+1)); % variance of d (exact)
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
    for i=1:np
        predInt(i,1)=residStdErrOLS^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorOLS(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(np)),': ',num2str(bestFitUnEq(np,1)),' deg; 95% prediction interval: ',num2str(bestFitUnEq(np,1)-predInt(np,1)),', ',num2str(bestFitUnEq(np,1)+predInt(np,1))]);
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
        for i=1:np
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n/2
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


disp('-------------------------------------------------------------------------');

if usingOctave
    fflush(stdout);
end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First find kd = rd/ro where rd = effective distance %
    % between the datasets, equation (4.7).               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=0;
    sumSameTrds=0;
    for i=1:n-1
        for j=i+1:n
            if abs(x(i,2)-x(j,2))<0.001
                % same time, different dataset
                m=m+1;
                sumSameTrds=sumSameTrds+(residualsOLS(i,1)-residualsOLS(j,1))^2;
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
    DF=n-k-1;
    disp(['kd=rd/ro = ',num2str(kd),', ',num2str(m),' co-located pairs->',num2str(n-k-1),'-',num2str(round(m*exp(-kd))),' = ',num2str(DF),' degrees of freedom']);
    
    if usingOctave
        fflush(stdout);
    end



disp('-------------------------------------------------------------------------');

if usingOctave
    fflush(stdout);
end

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tanh adjusted Durbin-Watson %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ed and Vd for un-transformed data already calculated
    atanhRho=(2/(nr-k1-3))*sqrt((nr-k1+3)/Vd)*(atanh(1-d/2)-atanh(1-Ed/2));
    rhoADW=real(tanh(atanhRho));

    roDW=-effDist/log(rhoDW);
    disp(['Durbin-Watson rho = ',num2str(rhoDW),' -> ro = -medianDist/ln(rho) = ',num2str(roDW)]);
    roADW=-effDist/log(rhoADW);
    disp(['tanh adjusted DW rho = ',num2str(rhoADW),' -> ro = -medianDist/ln(rho) = ',num2str(roADW)]);
    
disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=roDW;
    rd=kd*ro;
    disp(['GLS using ro = ',num2str(ro),' km, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            % 1D distance on time-line
            r=abs(x(i,2)-x(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
                %phi=exp(-(r/ro)-kd);
            end
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    %disp(['Inverting S']);
    if usingOctave
        fflush(stdout);
    end
    Sinv=inv(S);

    invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            % 1D distance on time-line (reduced dataset)
            r=abs(rx(i,2)-rx(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
            %    phi=exp(-(r/ro)-kd);
            end
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    invrXdSirX=inv(transpose(rx)*SRinv*rx);
    
betaGLSdw=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSdw;
residStdErrDW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
stdErrorDW=0;
tValue=0;
for m=1:k1
    stdErrorDW(m,1)=residStdErrDW*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLSdw(m,1)/stdErrorDW(m,1);
end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k+1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLSdw(m,1),stdErrorDW(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrDW),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

% Durbin-Watson equivalent calculation on equalised residuals
% using P=Sinv^0.5 i.e. P is symmetrical matrix satisfying P'P = Sinv
P=Sinv^0.5;
residualsEq=P*residualsGLS;
% reduce to nr residuals by premultiplying by R
residualsEq=R*residualsEq;
SSR=transpose(residualsEq)*residualsEq;
SSRFD=transpose(residualsEq)*A*residualsEq;
d=SSRFD/SSR;

% work out confidence limits of d
Pr=real(SRinv^0.5); % is real, but infintesimal imaginary component can cause problems
p=trace(A)-trace(transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX);
Ed=p/(nr-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(rx)*transpose(Pr)*A*A*Pr*rx*invrXdSirX)+trace((transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX)^2);
Vd=2*(q-p*Ed)/((nr-k-1)*(nr-k+1)); % variance of d (exact)
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
    tValue1=tinv(0.975,DF);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),DF);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:np
        predInt(i,1)=residStdErrDW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorDW(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(np)),': ',num2str(bestFitDW(np,1)),' deg; 95% prediction interval: ',num2str(bestFitDW(np,1)-predInt(np,1)),', ',num2str(bestFitDW(np,1)+predInt(np,1))]);
    end

    if plotBestFit==2
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
        for i=1:np
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n/2
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


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLS calculation by matrix algebra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=roADW;
    rd=kd*ro;
    disp(['GLS using ro = ',num2str(ro),' km, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            % 1D distance on time-line
            r=abs(x(i,2)-x(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
                %phi=exp(-(r/ro)-kd);
            end
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    %disp(['Inverting S']);
    if usingOctave
        fflush(stdout);
    end
    Sinv=inv(S);

invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            % 1D distance on time-line (reduced dataset)
            r=abs(rx(i,2)-rx(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
            %    phi=exp(-(r/ro)-kd);
            end
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    invrXdSirX=inv(transpose(rx)*SRinv*rx);
    
betaGLStadw=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLStadw;
residStdErrTADW=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
stdErrorTADW=0;
tValue=0;
for m=1:k1
    stdErrorTADW(m,1)=residStdErrTADW*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLStadw(m,1)/stdErrorTADW(m,1);
end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k+1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLStadw(m,1),stdErrorTADW(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrTADW),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

% Durbin-Watson equivalent calculation on equalised residuals
% using P=Sinv^0.5 i.e. P is symmetrical matrix satisfying P'P = Sinv
P=Sinv^0.5;
residualsEq=P*residualsGLS;
% reduce to nr residuals by premultiplying by R
residualsEq=R*residualsEq;
SSR=transpose(residualsEq)*residualsEq;
SSRFD=transpose(residualsEq)*A*residualsEq;
d=SSRFD/SSR;

% work out confidence limits of d
Pr=real(SRinv^0.5); % is real, but infintesimal imaginary component can cause problems
p=trace(A)-trace(transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX);
Ed=p/(nr-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(rx)*transpose(Pr)*A*A*Pr*rx*invrXdSirX)+trace((transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX)^2);
Vd=2*(q-p*Ed)/((nr-k-1)*(nr-k+1)); % variance of d (exact)
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
    tValue1=tinv(0.975,DF);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),DF);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:np
        predInt(i,1)=residStdErrTADW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorTADW(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(np)),': ',num2str(bestFitTADW(np,1)),' deg; 95% prediction interval: ',num2str(bestFitTADW(np,1)-predInt(np,1)),', ',num2str(bestFitTADW(np,1)+predInt(np,1))]);
    end

    if plotBestFit==3
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
        for i=1:np
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n/2
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
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k+1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLStadw(m,1),stdErrorExtrap(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end

    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:np
        predInt(i,1)=residStdErrTADW^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorExtrap(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(np)),': ',num2str(bestFitTADW(np,1)),' deg; 95% prediction interval: ',num2str(bestFitTADW(np,1)-predInt(np,1)),', ',num2str(bestFitTADW(np,1)+predInt(np,1))]);
    end

    if plotBestFit==4
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
        for i=1:np
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n/2
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

    disp('-------------------------------------------------------------------------');

    if usingOctave
        fflush(stdout);    
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Semivariogram exponential model fit %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot residual semivariogram
    gamma=0;
    h=0;
    m=0;
    hClass=0;
    gClass=0;
    nClass=0;
    i=0;
    for i=1:n-1
        for j=i+1:n
            m=m+1;
            gamma(m)=((residualsOLS(i)-residualsOLS(j))^2)/2;
            h(m)=abs(x(i,2)-x(j,2));
            %h(m)=sqrt((t(i)-t(j))^2+(u(i)-u(j))^2+(v(i)-v(j))^2); % 3D
            classNumb=round(h(m)+0.5);
            [~,classes]=size(nClass);
            if classNumb>classes
                nClass(classNumb)=1;
                gClass(classNumb)=gamma(m);
            else
                nClass(classNumb)=nClass(classNumb)+1;
                gClass(classNumb)=gClass(classNumb)+gamma(m);
            end
            hClass(classNumb)=classNumb-0.5;
        end
    end
    [~,classes]=size(nClass);
    for i=1:classes
        if nClass(i)>0;
            gClass(i)=gClass(i)/nClass(i);
        else
            gClass(i)=NaN;
        end
    end
    
    % weighted fit of exponential (residStdErrOLS^2)*[1-exp(-h/rs)] to semivariogram
    rs=0.1; % don't start too small, or improvement won't show up and loop will stop!
    residVar=residStdErrOLS^2;
    minSqdErr=0;
    for i=1:classes % normally for i=1:classes but here class 1 is h=0
        if nClass(i)>0
            %minSqdErr=minSqdErr+nClass(i)*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i) weighted
            minSqdErr=minSqdErr+(nClass(i)/hClass(i))*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i)/hClass(i) weighted
            %minSqdErr=minSqdErr+(nClass(i)/(hClass(i)^2))*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i)/(hClass(i)^2) weighted
        end
    end
    thisSqdErr=minSqdErr;
    lastSqdErr=minSqdErr*2;
    while minSqdErr<lastSqdErr
        rs=rs+0.01;
        lastSqdErr=thisSqdErr;
        thisSqdErr=0;
        for i=2:classes % normally for i=1:classes but here class 1 is h=0
            if nClass(i)>0;
                %thisSqdErr=thisSqdErr+nClass(i)*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i) weighted
                thisSqdErr=thisSqdErr+(nClass(i)/hClass(i))*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i)/hClass(i) weighted
                %thisSqdErr=thisSqdErr+(nClass(i)/(hClass(i)^2))*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2); % nClass(i)/(hClass(i)^2) weighted
            end
        end
        if thisSqdErr<=minSqdErr
            minSqdErr=thisSqdErr;
        end
    end
    
    % sigma^2 line
    dist(1)=0;
    dist(2)=max(h);
    rse2(1)=residStdErrOLS^2;
    rse2(2)=rse2(1);
    
    % variance*[1-exp(-h/rs)]
    hCurve=0;
    gCurve=0;
    for i=1:round(dist(2))
        hCurve(i)=(i-1);
        gCurve(i)=rse2(1)*(1-exp(-hCurve(i)/rs));
    end
    
%    if n<100
%        plot(h,gamma,'g.',hClass,gClass,'b-',dist,rse2,'r-',hCurve,gCurve,'m-');
%        legend('\gamma(i,j) points','\gamma in 1 km classes','\sigma^2',['semivariogram fit rs=',num2str(rs),'km']);
%    else
%        plot(hClass,gClass,'b-',dist,rse2,'r-',hCurve,gCurve,'m-');
%        legend('\gamma in 1 km classes','\sigma^2',['semivariogram fit ro=',num2str(rs),'km']);
%    end
%    xlabel('distance, km');
%    ylabel('\gamma(i,j), dB^2');
%    title(['OLS residual semivariogram, ',num2str(n),' points, ',num2str(k),' variables']);
%    disp(['semivariogram fitted ro = ',num2str(rs)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra using semivariogram fit %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=rs;
    rd=kd*ro;
    disp(['Semivariogram ro = ',num2str(ro),' km, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            % 1D distance on time-line
            r=abs(x(i,2)-x(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
                %phi=exp(-(r/ro)-kd);
            end
            S(i,j)=phi;
            S(j,i)=S(i,j);
        end
    end
    %disp(['Inverting S']);
    if usingOctave
        fflush(stdout);
    end
    Sinv=inv(S);

    invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            % 1D distance on time-line (reduced dataset)
            r=abs(rx(i,2)-rx(j,2));
            if abs(ds(i)-ds(j))<0.5
                % same dataset
                phi=exp(-r/ro);
            else
                % different datasets
                phi=exp(-sqrt((r/ro)^2+kd*kd));
            %    phi=exp(-(r/ro)-kd);
            end
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    invrXdSirX=inv(transpose(rx)*SRinv*rx);
    
betaGLSsvg=invXdSiX*transpose(x)*Sinv*y;
residualsGLS=y-x*betaGLSsvg;
residStdErrSVG=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
stdErrorSVG=0;
tValue=0;
for m=1:k1
    stdErrorSVG(m,1)=residStdErrSVG*sqrt(invXdSiX(m,m));
    tValue(m,1)=betaGLSdw(m,1)/stdErrorSVG(m,1);
end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k+1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLSsvg(m,1),stdErrorSVG(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrSVG),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

% Durbin-Watson equivalent calculation on equalised residuals
% using P=Sinv^0.5 i.e. P is symmetrical matrix satisfying P'P = Sinv
P=Sinv^0.5;
residualsEq=P*residualsGLS;
% reduce to nr residuals by premultiplying by R
residualsEq=R*residualsEq;
SSR=transpose(residualsEq)*residualsEq;
SSRFD=transpose(residualsEq)*A*residualsEq;
d=SSRFD/SSR;

% work out confidence limits of d
Pr=real(SRinv^0.5); % is real, but infintesimal imaginary component can cause problems
p=trace(A)-trace(transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX);
Ed=p/(nr-k-1); % mean of d (exact)
q=trace(A*A)-2*trace(transpose(rx)*transpose(Pr)*A*A*Pr*rx*invrXdSirX)+trace((transpose(rx)*transpose(Pr)*A*Pr*rx*invrXdSirX)^2);
Vd=2*(q-p*Ed)/((nr-k-1)*(nr-k+1)); % variance of d (exact)
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

    bestFitSVG=xp*betaGLSsvg;
    tValue1=tinv(0.975,DF);
    tValue2=tinv(1-0.025*(1-2*abs(betacdf(d/4,a,b)-0.5)),DF);
        % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:np
        predInt(i,1)=residStdErrSVG^2;
        for m=1:k1
            predInt(i,1)=predInt(i,1)+(stdErrorSVG(m,1)*dxp(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end
    if predictYears>0
        disp(['prediction for ',num2str(tp(np)),': ',num2str(bestFitSVG(np,1)),' deg; 95% prediction interval: ',num2str(bestFitSVG(np,1)-predInt(np,1)),', ',num2str(bestFitSVG(np,1)+predInt(np,1))]);
    end

    if plotBestFit==5
        outsideIP=0;
        for i=1:n/2
            if y1(i,1)>bestFitSVG(i,1)+predInt(i,1) || y1(i,1)<bestFitSVG(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
            if y2(i,1)>bestFitSVG(i,1)+predInt(i,1) || y2(i,1)<bestFitSVG(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y1,'b.-',t,y2,'k.-',tp,bestFitSVG,'r-',tp,bestFitSVG+predInt,'g-',tp,bestFitSVG-predInt,'g-');
        legend('NOAAGlobalTemp','HADCrut4','best fit (DW)','prediction 97.5%','prediction 2.5%','Location','northwest');
        %title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', good: Pr(>|d-',num2str(Ed,3),'|): ',num2str(1-2*abs(betacdf(d/4,a,b)-0.5)),']']);
        title([num2str(n),' measured data points, and best fit, model order = ',num2str(modelOrder),', ',dataName,' data  [after GLS transform, Durbin-Watson: ',num2str(d),', [',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel('Year');
        ylabel('relative global temperature, land & sea, deg C');
        % output to csv
        CSVoutput=0;
        for i=1:np
            CSVoutput(i,1)=tp(i);                         % year AD
            if i>n/2
                CSVoutput(i,2)=NaN;
                CSVoutput(i,3)=NaN;
            else
                CSVoutput(i,2)=y1(i,1);                    % observed data
                CSVoutput(i,3)=y2(i,1);
            end
            CSVoutput(i,4)=bestFitSVG(i,1);               % predicted value
            CSVoutput(i,5)=bestFitSVG(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,6)=bestFitSVG(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'Order',num2str(modelOrder),'polynomialTrendSVG.csv'],CSVoutput);
    end
        
end


disp('=========================================================================');

