T = 5;
d = 2;
numpoints = 50;
betaActual = zeros(1,T)
FK= zeros(1,T)
betaTheory = zeros(1,T);
betaTheory_sample = zeros(1,T);
betaTheory_sample2 = zeros(1,T);
percentsep = zeros(1,T);
ahist = zeros(2*numpoints,d,T);
bhist = zeros(2*numpoints,d,T);
a1 = zeros(1,T);
a2= zeros(1,T);
a3=zeros(1,T);
a4= zeros(1,T);
linem = zeros(1,T);
lineb= zeros(1,T);
plotsingle = 0;
gen = 1;
for d = 2:2
for trial = 1:T
if gen == 1 
    
[a,b,linem,lineb] = dividePoints(numpoints,d);
%scatter(a(:,1),a(:,2)); hold on; 
%scatter(b(:,1),b(:,2)); hold on;
%x = linspace(-25,25);
%y = linem*x+lineb;
%plot(x,y);
for indx = 1:T
    atemp = a(indx,:);
    btemp = b(indx,:);
    a(indx,:) = btemp;
    b(indx,:) = atemp;
end

A = size(a,1);
B = size(b,1);

ahist(1:A,:,trial) = a;
bhist(1:B,:,trial) = b;
end


X = [a;b];
Y= [ zeros(A,1);ones(B,1)];
U = [X Y];
[betaActual(1,trial), linem(1,trial), lineb(1,trial)] = weakseparator(a,b,d);
%[betaActual(1,trial), a1(1,trial),a2(1,trial),a3(1,trial),a4(1,trial)] = weakseparator(a,b,d);
betaActual(1,trial) = betaActual(1,trial)/(A+B);


[C,Cn] = dptuples(U,d);
comblength = size(C,3);
separableTuples = 0;
sampleSize= round((.1+.1*t)*comblength); 

%%%%%Setup for Sampling Test Results%%%%%%%%
separability = zeros(1,comblength);
sampleseparability = zeros(1,sampleSize);
samples= datasample(1:comblength,sampleSize,'Replace',false);
separableSamples = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Setup for Sampling point subsets%%%%%%%%%
sampleSize2 = round((.1+.1*t)*(A+B));
samplePoints2 = datasample(1:(A+B),sampleSize2,'Replace',false);
sampleTuples2 = combnk(samplePoints,d+1);
separability2 = zeros(1,comblength);
sampleseparability2 = zeros(1,size(sampleTuples2,1));
Cn = sort(Cn,2);
sampleTuples2 = sort(sampleTuples,2);
sampleTuples2 = sortrows(sampleTuples);
separableSamples2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for dir = -1:2:1
    j=1;
    j2 =1;
for i = 1:comblength
    aindex= find(C(:,d+1,i) == 0);
    bindex = find(C(:,d+1,i) == 1);
    if (aindex)
        atest = C(aindex,1:d,i)
    else
        atest = [];
    end
    if(bindex)
        btest = C(bindex,1:d,i)  
    else
        btest = [];
    end
    separability(1,i) = testsep(atest,btest,d,dir);
    if any(samples(1,:) == i)
        sampleseparability(1,j) = separability(1,i);
        j = j+1;
    end
    if ismember(Cn(i,:),sampleTuples2,'rows','legacy')
        sampleseparability2(1,j2) = separability(1,i);
        j2 = j2+1;
    end
end
if (sum(separability(1,:))>separableTuples)
    separableTuples = sum(separability(1,:));
end
if (sum(sampleseparability(1,:)) > separableSamples)
    separableSamples = sum(sampleseparability(1,:));
end
if (sum(sampleseparability2(1,:)) > separableSamples2)
    separableSamples2 = sum(sampleseparability2(1,:));
end
end
percentsep(1,trial) = separableTuples/comblength;
alpha = separableTuples/comblength;
alpha_sample = separableSamples/sampleSize;
alpha_sample2 = separableSamples2/size(sampleTuples2,1);
FK(1,trial) = alpha;
betaTheory(1,trial) = (rFinder(separableTuples,(A+B),d)+d+1)/(A+B);
betaTheory_sample(1,trial) = (rFinder(alpha_sample*comblength,(A+B),d)+d+1)/(A+B);
betaTheory_sample2(1,trial) = (rFinder(separableSamples,size(samplePoints,2),d)+d+1)/(size(samplePoints2,2));

save('FK_Sampling_'+string(A+B)+'_T_'+string(T));
if plotsingle == 1;
    scatter(ahist(:,1,T),ahist(:,2,T)); hold on; scatter(bhist(:,1,T),bhist(:,2,T));hold on;
    x = linspace(-20,20);
    y = linem(T)*x+lineb(T);
    plot(x,y);
    %%%3dplotting
    plot3(ahist(:,1,T),ahist(:,2,T),ahist(:,3,T),'o'); hold on; plot3(bhist(:,1,T),bhist(:,2,T),bhist(:,3,T),'x');
    [x y] = meshgrid(-20:1:20);
    hold on;
    z = -a1(T)/a3(T).*x-a2(T)/a3(T).*y+a4(T)/a3(T);
    surf(x,y,z);
end

%[r pk fk] = rFinder(alpha,(A+B),d);
end
end


function [C] = dpptuples(U,d)
length = size(U,1);    
Cn = combnk(1:length,d+2);
Cnsize = size(Cn,1);
C = zeros(4,size(U,2),Cnsize)
for i = 1:Cnsize
    C(1,:,i) = U(Cn(i,1),:);
    C(2,:,i) = U(Cn(i,2),:);
    C(3,:,i) = U(Cn(i,3),:);
    C(4,:,i) = U(Cn(i,4),:);
end
end
function[C] = dptuples(U,d,p)
length = size(U,1);
Cn = combnk(1:length,d+1);
Cnsize = size(Cn,1);
C = zeros(d+1,size(U,2),Cnsize);
for i = 1:Cnsize
    for j = 1:d+1
    C(j,:,i) = U(Cn(i,j),:);
    end
   
end
end
function [separable] = testsep(a,b,d,dir)
sizea = size(a,1);
sizeb = size(b,1);
if sizea== d+1 || sizeb == d+1
    separable = 1;
   
else
A1 = [-a; b]
bs1 = [ones(1,sizea) -ones(1,sizeb)];
A2= [a;-b]
bs2 = [-ones(1,sizea) ones(1,sizeb)];
f = zeros(1,d);
if dir >0
separable = any(logical(linprog(f,A1,bs1)));
else
separable = any(logical(linprog(f,A2,bs2)));
end
%separable = separable1||separable2;
end
end

function [a,b] = uniform(A,B,d)
a = zeros(A,d);
b = zeros(B,d);
for n = 1:A
    r = zeros(1,d);
   for i = 1:d
       r(1,i) = ((rand()-.5)*50);
   end
a(n,:) =r;
end
for n = 1:B
r = zeros(1,d);
    for i = 1:d
        r(1,i) = ((rand()-.5)*50);
    end
b(n,:) =r;
end    
end

function [a,b]= normal(A,B,d,Amean,Bmean,Astd,Bstd)

a = zeros(A,d);
b = zeros(B,d);
for n = 1:A
    a(n,:) = Astd*randn(1,d)+Amean;
end
for n = 1:B
    b(n,:) = Bstd*randn(1,d)+Bmean;
end    
end

function [a,b] = setBeta(A,B,beta,offset)
n = A+B;
v = floor(beta*n/2);
q = 3;
for i = 1:v
a(i,:) = [-3*cos(2*pi/v*i)+10 3*sin(2*pi/v*i)];
b(i,:) = [3*cos(2*pi/v*i)+20 3*sin(2*pi/v*i)];
end
for i = v+1:A
a(i,:) = [cos(2*pi/(A-v)*i)+20 sin(2*pi/(A-v)*i)-offset];
b(i,:) = [cos(2*pi/(B-v)*i)+10 sin(2*pi/(B-v)*i)-offset];
a(i,:) = [18 offset];
b(i,:) = [12 offset];
%q = 6;
%a(i,:) = [(.5*b(q,1)+.5*b(q+1,1))+.01 (.5*b(q,2)+.5*b(q+1,2))-.01];
%q = 3;
%b(i,:)= [(.65*a(q,1)+.35*a(q+1,1))-.01 (.65*a(q,2)+.35*a(q+1,2))-.01];

end
end
function [a,b] = setNoisyBeta(A,B,beta)
n = A+B;
v = floor(beta*n/2);

for i = 1:v
a(i,:) = [3*cos(2*pi/v*i)-10+(rand()-.5) 3*sin(2*pi/v*i)+(rand()-.5)];
b(i,:) = [3*cos(2*pi/v*i)+10+(rand()-.5) 3*sin(2*pi/v*i)+(rand()-.5)];
end
for i = v+1:A
a(i,:) = [cos(2*pi/(A-v)*i)+10+(rand()-.5) sin(2*pi/(A-v)*i)+(rand()-.5)];
b(i,:) = [cos(2*pi/(B-v)*i)-10+(rand()-.5) sin(2*pi/(B-v)*i)+(rand()-.5)];
end
end
function rmax = rFinder(alpha,n,d)
k = d+1;
%limit = alpha*nchoosek(n,k);
limit = alpha;
for r = k:n-d
    pk = 0;
    for i = 0:d
        pk = pk+nchoosek(r,(k-i))*nchoosek(n-r,i);
    end
    if limit <= pk
        rmax = r-1;
        break;
    end
end
end
function [a,b] = randpos_labels(numpoints,d,p)

a = [];
b = [];
for n = 1:numpoints
    r = rand();
    t = zeros(1,d);
   for i = 1:d
       t(1,i) = ((rand()-.5)*50);
   end
   if r< p
       a = [a ;t];
   else
       b = [b ;t];
    end
  
end
end
function [a,b,linem,lineb] = dividePoints(n,d)
linem = randi(10)/randi(10);
lineb = randi(10)-6;
a = zeros(n,2);
b = zeros(n,2);
for i = 1:n
    a(i,1) = (rand()-.5)*50;
    b(i,1) = (rand()-.5)*50;
    a(i,2) = linem*a(i,1)+lineb+rand()*(50-linem*a(i,1)-lineb);
    b(i,2) = linem*b(i,1)+lineb-rand()*(50-linem*b(i,1)-lineb);
end
end
