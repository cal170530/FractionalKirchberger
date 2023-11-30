T = 5000;
betaTheory = betaTheory(1,1:T);
betaActual= betaActual(1,1:T);
FK = FK(1,1:T);

xn = size(betaActual,2);
betaActualsorted = zeros(size(betaActual));
betaTheorysorted = zeros(size(betaTheory));
alphasorted = sort(FK);
alphaxlabels = [];
alphaunique = unique(alphasorted);
binning = 0;
numBins = 20;
%{
for i = 1:size(alphasorted,2)
    if mod(i,100) == 0
        alphaxlabels(1,i/100) = alphasorted(1,i);
    end
    for j = 1:size(alphasorted,2)
        if alphasorted(1,i) == FK(1,j)
            betaActualsorted(1,i) = betaActual(1,j);
            betaTheorysorted(1,i) = betaTheory(1,j);
        end
    end
end
%}
for i = 1:size(alphaunique,2)
 if mod(i,10 == 0)
     alphaxlabels = [alphaxlabels string(alphaunique(1,i))];
 end
end
betaActualAvg = zeros(1,size(alphaunique,2));
betaTheoryAvg = zeros(1,size(alphaunique,2));

for i = 1:size(alphaunique,2)
 numinstances = sum(FK(1,:)== alphaunique(1,i));
 betaTheoryAvg(1,i) = sum(betaTheory.*(FK(1,:) == alphaunique(1,i)))/numinstances;
 betaActualAvg(1,i) = sum(betaActual.*(FK(1,:) == alphaunique(1,i)))/numinstances;
end

if binning ==1
    Y = discretize(alphaunique, numBins);
    betaTheoryBinned = zeros(1,numBins);
    betaActualBinned = zeros(1,numBins);
    alphaUniqueBinned = zeros(1,numBins);
    for n = 1:numBins
        num_bin_occurances = sum(Y(1,:) == n);
        alphaUniqueBinned(1,n) = sum(alphaunique(1,:).*(Y(1,:) == n))/num_bin_occurances;
        betaTheoryBinned(1,n) = sum(betaTheoryAvg(1,:).*(Y(1,:) == n))/num_bin_occurances;
        betaActualBinned(1,n) = sum(betaActualAvg(1,:).*(Y(1,:) == n))/num_bin_occurances;
    end
    alphaunique = alphaUniqueBinned;
    betaTheoryAvg = betaTheoryBinned;
    betaActualAvg = betaActualBinned;
end
figure;
f = figure;
f.Position=[10 10 1000 800];
%data = [betaActualsorted' betaTheorysorted'];
data = [betaTheoryAvg' betaActualAvg'];
%x= 1:xn;
x = 1:size(alphaunique,2);
plot(alphaunique,betaActualAvg,'b',alphaunique,betaTheoryAvg,'r');
avgerror = sum(betaActual-betaTheory)/xn;

%bar(alphaunique,data);
%plot(x,betaActualBinned,'b',x,betaTheoryBinned,'r');
%plot(alphasorted,betaActualsorted,'b',alphasorted,betaTheorysorted,'r');
axis([-inf inf 0 1] );
title('d= 3, A = 10, B = 10 (Uniform)');
xlabel('\alpha');
ylabel('\beta (avg error = '+string(avgerror)+')');
xticks(.75:.0125:1);
%set(gca,'xticklabel',alphaxlabels)
legend('Actual','FKLowerBound');
figure;
hist(FK,size(alphaunique,2));
xlabel('\alpha');
ylabel('number of occurences');

alphaplot = alphaunique(1,1:3:end);
betaplot = betaActualAvg(1,1:3:end);
for i = 1:size(alphaplot,2)
    p = alphaplot(1,i);
    p1 = p-3*sqrt(10*p*(1-p))