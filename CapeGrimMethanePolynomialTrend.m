%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads Cape Grim Methane data, then does OLS and GLS linear regression. %
% GLS assumes residual correlation is by exponential functions of distance.           % 
% This code is for Matlab or Octave (clear or set the usingOctave flag).              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usingOctave=false;  % forces flushing of the display buffer in Octave to see progress (not for Matlab)

modelOrder=3;     % regression model order (max 9 with this code)

% Read data:
methane=csvread('CapeGrimMethane.csv');

plotBestFit=0; % plot best fit with PI: 0=no, 1=OLS, 2=DW, 3=TADW, 4=DW-TADW extrapolated, 5=semivariogram
% if 0, semivariogram will be plotted instead of data with best fir and PI

doSemivar=true; % do GLS based on semivariogram of OLS residuals

[n,~]=size(methane);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy the first 5 records as is, then take means 12 at a time of the rest %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=zeros(5,n);
for i=1:5
    M(i,i)=1;
end
for i=6:n
    j=round((i-6-mod(i-6,12))/12)+6;
    M(j,i)=1/12;
end
lastRowSum=sum(M(j,:)); % make the last row sum = 1
M(j,:)=M(j,:)/lastRowSum;
methaneA=M*methane;
[n,~]=size(methaneA);

%%%%%%%%%%%%%%%%%%%%%
% Assemble the data %
%%%%%%%%%%%%%%%%%%%%%
dataMeans=mean(methaneA);
meanT=dataMeans(1);
m=0;
t=0;  % time
x=0;  % OLS/GLS matrix
y=0;  % correlated data
for i=1:n
    dateDecimal=methaneA(i,1);
    CH4ppb=methaneA(i,2);
    t(i,1)=dateDecimal;
    x(i,1)=1;
    %ds(i)=1; % all the one dataset
    if modelOrder>0
        x(i,2)=t(i)-meanT;
    end
    if modelOrder>1
        for j=2:modelOrder
            x(i,j+1)=x(i,2)^j;
        end
    end
    y(i,1)=CH4ppb;
end

% decorrelate the explanatory variables for accurate prediction intervals
if modelOrder>1
    for k=2:modelOrder
        % remove any correlation between x(:,k+1) and x(:,1:k) by subtracting
        % an OLS estimate of x(:,k+1) in terms of x(:,1:k)
        X=x(:,1:k);
        Y=x(:,k+1);
        betaDecorr=inv((transpose(X)*X))*transpose(X)*Y;
        x(:,k+1)=Y-X*betaDecorr;
    end
end

[n,k1]=size(x); % k1 = k + 1 = no. of variables + intercept
k=k1-1;
% calculate dx = Xi - Xmean
meanX=mean(x);
meanX(1,1)=0; % except dx(intercept) = 1
dx=meanX*0;
for i=1:n
    dx(i,:)=x(i,:)-meanX;
end
dataName=['CapeGrimCH4GLS',num2str(n),'points'];


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

disp('=========================================================================');
disp('Data: Atmospheric background methane at Cape Grim');
disp('from http://www.csiro.au/greenhouse-gases/');
if modelOrder>1
    disp('The 2nd and higher order expanatory variables are decorrelated w.r.t. previous variables');
end
disp('-------------------------------------------------------------------------');
disp('OLS calculation by matrix algebra:');
invXdX=inv(transpose(x)*x);
betaOLS=invXdX*transpose(x)*y;
residualsOLS=y-x*betaOLS;
residStdErrOLS=sqrt(transpose(residualsOLS)*residualsOLS/(n-k1));
stdErrorOLS=0;
tValue=0;
for m=1:k1
    stdErrorOLS(m,1)=residStdErrOLS*sqrt(invXdX(m,m));
    tValue(m,1)=betaOLS(m,1)/stdErrorOLS(m,1);
end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaOLS(m,1),stdErrorOLS(m,1),tValue(m,1),n-k1,2*(1-tcdf(abs(tValue(m,1)),n-k-1)));
end
rmsError=sqrt((transpose(residualsOLS)*residualsOLS)/n);
disp(['Residual standard error: ',num2str(residStdErrOLS),' on ',num2str(n-k1),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

    

if doSemivar
    % plot residual semivariogram
    gamma=0;
    h=0;
    m=0;
    hClass=0;
    gClass=0;
    nClass=0;
    i=0;
    for i=1:n
        for j=i:n
            m=m+1;
            gamma(m)=((residualsOLS(i)-residualsOLS(j))^2)/2;
            h(m)=abs(t(i)-t(j));
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
        if nClass(i)>0
            gClass(i)=gClass(i)/nClass(i);
        else
            gClass(i)=NaN;
        end
    end
    % weighted fit of exponential (residStdErrOLS^2)*[1-exp(-h/rs)] to semivariogram
    rs=0.1;
    residVar=residStdErrOLS^2;
    minSqdErr=0;
    for i=1:classes
        if nClass(i)>0
            minSqdErr=minSqdErr+nClass(i)*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2);
        end
    end
    thisSqdErr=minSqdErr;
    lastSqdErr=minSqdErr*2;
    while minSqdErr<lastSqdErr
        rs=rs+0.01;
        lastSqdErr=thisSqdErr;
        thisSqdErr=0;
        for i=1:classes
            if nClass(i)>0
                thisSqdErr=thisSqdErr+nClass(i)*((gClass(i)-residVar*(1-exp(-hClass(i)/rs)))^2);
            end
        end
        if thisSqdErr<minSqdErr
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
    
    if plotBestFit==0
    
    if n<100
        plot(h,gamma,'g.',hClass,gClass,'b-',dist,rse2,'r-',hCurve,gCurve,'m-');
        legend('\gamma(i,j) points','\gamma in 100 km classes','\sigma^2',['semivariogram fit rs=',num2str(rs),'years']);
    else
        plot(hClass,gClass,'b-',dist,rse2,'r-',hCurve,gCurve,'m-');
        legend('\gamma in 100 km classes','\sigma^2',['semivariogram fit ro=',num2str(rs),'years']);
    end
    xlabel('distance, years');
    ylabel('\gamma(i,j), deg C^2');
    title(['OLS residual semivariogram, ',num2str(n),' samples, ',num2str(k-1),' variables']);
    
    end
    disp(['semivariogram fitted ro = ',num2str(rs),' years    Making list of neighbours ...']);
else
    disp('Making progressive list of nearest neighbours ...');
end
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make list of "adjacent" points for "forward differences" for generalised Durbin-Watson %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if usingOctave
    fflush(stdout);
end
% Find starting point with the greatest sum of distances to all other points
rSumMax=0;
nSumMax=0;
for i=1:n
    rSumThisOne=0;
    for j=1:n
        r=abs(t(i)-t(j));
        rSumThisOne=rSumThisOne+r;
    end
    if rSumThisOne>rSumMax
        rSumMax=rSumThisOne;
        nSumMax=i;
    end
end
%disp(['Starting point is no. ',num2str(nSumMax),' at lat ',num2str(lat(nSumMax)),', long ',num2str(long(nSumMax))]);
% also generate reduction matrix to pre-multiply x, y or e to average co-located points
R=zeros(1,n); % first row
Rrow=1;
nThisRow=1;
R(Rrow,nSumMax)=1;

adjacentPoint=zeros(1,n);
adjacentDist=zeros(1,n);
%adjacentDistToD=zeros(1,n);
thisPoint=nSumMax; % starting point
startPoint=nSumMax;
% find next point
rMin=999999999;
for j=1:n
    r=abs(t(thisPoint)-t(j));
    if j~=thisPoint && r<rMin
        % closest so far
        rMin=r;
        adjacentPoint(thisPoint)=j;
        adjacentDist(thisPoint)=r;
    end
end
% and then continue until no more nearest new neighbours found
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
            r=abs(t(thisPoint)-t(j));
            if j~=thisPoint && r<rMin
                % closest so far
                rMin=r;
                adjacentPoint(thisPoint)=j;
                adjacentDist(thisPoint)=r;
            end
        end
    end
end


% reduced coordinate list
rt=R*t;

% use premultiplication matrix R to merge co-loacted points and sort them
% in nearest new neighbour order
[nr,~]=size(R);
effDist=sum(adjacentDist)/(nr-1); % the last point on the path has adjacentDist = 0
disp(['Starting with no ',num2str(nSumMax),', time ',num2str(t(nSumMax)),', mean non-zero distance = ',num2str(effDist),' years']);
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
SSRFD=(transpose(R*residualsOLS)*A*(R*residualsOLS));
SSR=transpose(R*residualsOLS)*(R*residualsOLS);
d=SSRFD/SSR;

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

    bestFitOLS=x*betaOLS;
    tValue1=tinv(0.975,n-k-1);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n
        predInt(i,1)=residStdErrOLS^2;
        for m=1:k+1
            predInt(i,1)=predInt(i,1)+(stdErrorOLS(m,1)*dx(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end

    if plotBestFit==1
        outsideIP=0;
        for i=1:n
            if y(i,1)>bestFitOLS(i,1)+predInt(i,1) || y(i,1)<bestFitOLS(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y,'b.',t,bestFitOLS,'r-',t,bestFitOLS+predInt,'g-',t,bestFitOLS-predInt,'g-');
        legend('Observed','best fit (OLS)','prediction 97.5%','prediction 2.5%','Location','northwest');
        title([num2str(n),' measured data points, and best fit, ',dataName,' data  [GLS transformed residual Durbin-Watson: ',num2str(d),',  ',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel(varName(2,:));
        ylabel('A 0.01%, dB');
        % output to csv
        CSVoutput=0;
        for i=1:n
            CSVoutput(i,1)=t(i);                          % parameter
            CSVoutput(i,2)=y(i,1);                        % observed data
            CSVoutput(i,3)=bestFitOLS(i,1);              % predicted value
            CSVoutput(i,4)=bestFitOLS(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,5)=bestFitOLS(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'order',num2str(modelOrder),'polynomialOLS.csv'],CSVoutput);
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
            if abs(t(i)-t(j))<0.001
                % same location, different data point
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
        phi=0;
    end
    DF=n-k-1;
    disp(['kd = rd/ro = ',num2str(kd),' = (same location effective distance)/ro - see (4.7)']);
    
    if usingOctave
        fflush(stdout);
    end


disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Simple Durbin-Watson %
    %%%%%%%%%%%%%%%%%%%%%%%%
    rhoDW=1-d/2;
    % and tanh adjusted
    atanhRho=(2/(nr-k1-3))*sqrt((nr-k1+3)/Vd)*(atanh(1-d/2)-atanh(1-Ed/2));
    rhoADW=real(tanh(atanhRho));
if rhoDW>0
    roDW=-effDist/log(rhoDW);
    disp(['DW rho = 1 - d/2 = ',num2str(rhoDW),' -> ro = -effDist/ln(rho) = ',num2str(roDW)]);
    
disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=roDW;
    rd=kd*ro;
    disp(['GLS using ro = ',num2str(ro),' years, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            r=abs(t(i)-t(j));
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
    invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            r=abs(rt(i)-rt(j));
            phi=exp(-sqrt((r/ro)^2+kd*kd));
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    invrXdSirX=inv(transpose(rx)*SRinv*rx);
    
    betaGLS=invXdSiX*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLS;
    residStdErrGLS=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
    stdErrorGLS=0;
    tValue=0;
    for m=1:k1
        stdErrorGLS(m,1)=residStdErrGLS*sqrt(invXdSiX(m,m));
        tValue(m,1)=betaGLS(m,1)/stdErrorGLS(m,1);
    end
    stdErrorDW=stdErrorGLS; % } save for t
    betaGLSdw=betaGLS;      % } extrapolation
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLS(m,1),stdErrorGLS(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrGLS),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);


% Durbin-Watson equivalent calculation on transformed residuals
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

    bestFitGLS=x*betaGLS;
    tValue1=tinv(0.975,DF);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n
        predInt(i,1)=residStdErrGLS^2;
        for m=1:k+1
            predInt(i,1)=predInt(i,1)+(stdErrorGLS(m,1)*dx(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end

    if plotBestFit==2
        outsideIP=0;
        for i=1:n
            if y(i,1)>bestFitGLS(i,1)+predInt(i,1) || y(i,1)<bestFitGLS(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y,'b.',t,bestFitGLS,'r-',t,bestFitGLS+predInt,'g-',t,bestFitGLS-predInt,'g-');
        legend('Observed','best fit GLS DW','prediction 97.5%','prediction 2.5%','Location','northwest');
        title([num2str(n),' measured data points, and best fit, ',dataName,' data  [GLS transformed residual Durbin-Watson: ',num2str(d),',  ',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel(varName(2,:));
        ylabel('A 0.01%, dB');
        % output to csv
        CSVoutput=0;
        for i=1:n
            CSVoutput(i,1)=t(i);                          % parameter
            CSVoutput(i,2)=y(i,1);                        % observed data
            CSVoutput(i,3)=bestFitGLS(i,1);              % predicted value
            CSVoutput(i,4)=bestFitGLS(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,5)=bestFitGLS(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'order',num2str(modelOrder),'polynomialGLS-DW.csv'],CSVoutput);
    end

end


else
    disp(['simple DW rho = 1 - d/2 = ',num2str(rhoADW),' -> skip GLS estimation']);
end % if rhoDW>0


disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end


    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tanh adjusted Durbin-Watson %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rhoADW>0
    roADW=-effDist/log(rhoADW);
    disp(['tanh adjusted DW rho = ',num2str(rhoADW),' -> ro = -effDist/ln(rho) = ',num2str(roADW)]);
    
disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=roADW;
    rd=kd*ro;
    disp(['GLS using ro = ',num2str(ro),' years, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            r=abs(t(i)-t(j));
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
    invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            r=abs(rt(i)-rt(j));
            phi=exp(-sqrt((r/ro)^2+kd*kd));
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    rx=R*x;
    invrXdSirX=inv(transpose(rx)*SRinv*rx);
    
    betaGLS=invXdSiX*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLS;
    residStdErrGLS=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
    stdErrorGLS=0;
    tValue=0;
    for m=1:k1
        stdErrorGLS(m,1)=residStdErrGLS*sqrt(invXdSiX(m,m));
        tValue(m,1)=betaGLS(m,1)/stdErrorGLS(m,1);
    end
    stdErrorTADW=stdErrorGLS;
    betaGLStadw=betaGLS;
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLS(m,1),stdErrorGLS(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrGLS),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);


% Durbin-Watson equivalent calculation on transformed residuals
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

    bestFitGLS=x*betaGLS;
    tValue1=tinv(0.975,DF);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n
        predInt(i,1)=residStdErrGLS^2;
        for m=1:k+1
            predInt(i,1)=predInt(i,1)+(stdErrorGLS(m,1)*dx(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end

    if plotBestFit==3
        outsideIP=0;
        for i=1:n
            if y(i,1)>bestFitGLS(i,1)+predInt(i,1) || y(i,1)<bestFitGLS(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y,'b.',t,bestFitGLS,'r-',t,bestFitGLS+predInt,'g-',t,bestFitGLS-predInt,'g-');
        legend('Observed','best fit TADW','prediction 97.5%','prediction 2.5%','Location','northwest');
        title([num2str(n),' measured data points, and best fit, ',dataName,' data  [GLS transformed residual Durbin-Watson: ',num2str(d),',  ',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel(varName(2,:));
        ylabel('A 0.01%, dB');
        % output to csv
        CSVoutput=0;
        for i=1:n
            CSVoutput(i,1)=t(i);                          % parameter
            CSVoutput(i,2)=y(i,1);                        % observed data
            CSVoutput(i,3)=bestFitGLS(i,1);              % predicted value
            CSVoutput(i,4)=bestFitGLS(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,5)=bestFitGLS(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'order',num2str(modelOrder),'polynomialGLS-TADW.csv'],CSVoutput);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply t and Sb extrapolation since d test ok %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('------- Regression coeff std errors and t extrapolated: 2*TADW-DW -------');
    stdErrorExtrap=0;
    tValue=0;
    for m=1:k+1
        stdErrorExtrap(m,1)=2*stdErrorTADW(m,1)-stdErrorDW(m,1);
        tValue(m,1)=2*(betaGLStadw(m,1)/stdErrorTADW(m,1))-(betaGLSdw(m,1)/stdErrorDW(m,1));
    end
    disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
    for m=1:k1
        fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLStadw(m,1),stdErrorExtrap(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
    end

    tValue1=tinv(0.975,DF);
    % predInt is prediction interval based on regression coeff std errors
    predInt=0;
    for i=1:n
        predInt(i,1)=residStdErrGLS^2;
        for m=1:k+1
            predInt(i,1)=predInt(i,1)+(stdErrorExtrap(m,1)*dx(i,m))^2;
        end
        predInt(i,1)=tValue1*sqrt(predInt(i,1));
    end

    if plotBestFit==4
        outsideIP=0;
        for i=1:n
            if y(i,1)>bestFitGLS(i,1)+predInt(i,1) || y(i,1)<bestFitGLS(i,1)-predInt(i,1)
                outsideIP=outsideIP+1;
            end
        end
        plot(t,y,'b.',t,bestFitGLS,'r-',t,bestFitGLS+predInt,'g-',t,bestFitGLS-predInt,'g-');
        legend('Observed','best fit DW-TADW extrap','prediction 97.5%','prediction 2.5%','Location','northwest');
        title([num2str(n),' measured data points, and best fit, ',dataName,' data  [GLS transformed residual Durbin-Watson: ',num2str(d),',  ',num2str(100*outsideIP/n),'% ouside IP]']);
        xlabel(varName(2,:));
        ylabel('A 0.01%, dB');
        % output to csv
        CSVoutput=0;
        for i=1:n
            CSVoutput(i,1)=t(i);                          % parameter
            CSVoutput(i,2)=y(i,1);                        % observed data
            CSVoutput(i,3)=bestFitGLS(i,1);              % predicted value
            CSVoutput(i,4)=bestFitGLS(i,1)+predInt(i,1); % upper prediction interval
            CSVoutput(i,5)=bestFitGLS(i,1)-predInt(i,1); % lower prediction interval
        end
        csvwrite([dataName,'order',num2str(modelOrder),'polynomialGLS-DW-TADWextrap.csv'],CSVoutput);
    end


end


else
    disp(['tanh adjusted DW rho = ',num2str(rhoADW),' -> skip GLS estimation']);
end % if rhoADW>0
    

if doSemivar

disp('-------------------------------------------------------------------------');

    
    if usingOctave
        fflush(stdout);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLS calculation by matrix algebra using semivariogram fit %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro=rs;
    rd=kd*ro;
    disp(['Semivariogram ro = ',num2str(ro),' years, same location effective distance rd = ',num2str(rd),' km']);
    if usingOctave
        fflush(stdout);
    end
    % calculate S first, then invert (hopefully)
    S=eye(n);
    for i=1:n-1
        for j=i+1:n
            r=abs(t(i)-t(j));
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
    invXdSiX=inv(transpose(x)*Sinv*x);

    % repeat for reduced dataset
    Sr=eye(nr);
    for i=1:nr-1
        for j=i+1:nr
            r=abs(rt(i)-rt(j));
            phi=exp(-sqrt((r/ro)^2+kd*kd));
            Sr(i,j)=phi;
            Sr(j,i)=S(i,j);
        end
    end
    SRinv=inv(Sr);
    rx=R*x;
    invrXdSirX=inv(transpose(rx)*SRinv*rx);


    betaGLS=invXdSiX*transpose(x)*Sinv*y;
    residualsGLS=y-x*betaGLS;
    residStdErrGLS=sqrt(transpose(residualsGLS)*Sinv*residualsGLS/(DF));
    stdErrorGLS=0;
    tValue=0;
    for m=1:k1
        stdErrorGLS(m,1)=residStdErrGLS*sqrt(invXdSiX(m,m));
        tValue(m,1)=betaGLS(m,1)/stdErrorGLS(m,1);
    end
disp('                 Estimate     Std. Error    t value      DF      Pr(>|t|)');
for m=1:k1
    fprintf('%s  %12.6f  %12.6f  %9.3f  %7.0f  %12.6f\n',varName(m,:),betaGLS(m,1),stdErrorGLS(m,1),tValue(m,1),DF,2*(1-tcdf(abs(tValue(m,1)),DF)));
end
rmsError=sqrt((transpose(residualsGLS)*residualsGLS)/n);
disp(['Residual standard error: ',num2str(residStdErrGLS),' on ',num2str(DF),' Deg.Fr.  RMS error: ',num2str(rmsError)]);

% Durbin-Watson equivalent calculation on transformed residuals
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
end


end
    
disp('=========================================================================');
    
