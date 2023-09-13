T = 10000;
for numpoints = 10
for d = 3:3
%d = 2;
%numpoints = 10;
betaActual = zeros(1,T)
FK= zeros(1,T)
betaTheory = zeros(1,T);
percentsep = zeros(1,T);
ahist = zeros(numpoints,d,T);
bhist = zeros(numpoints,d,T);
a1 = zeros(1,T);
a2= zeros(1,T);
a3=zeros(1,T);
a4= zeros(1,T);
linem = zeros(1,T);
lineb= zeros(1,T);
plotsingle = 0;

%for d = 2:2
for trial = 1:T

A = numpoints;
B = numpoints;
%[a,b] = setBeta(A,B,.45+.05*trial);
%[a,b] = uniform(A,B,d);
[a,b] = randpos_labels(20,d,.5);
%ahist(:,:,trial) = a;
%bhist(:,:,trial) = b;
A = size(a,1);
B = size(b,1);
X = [a;b];
Y= [ zeros(A,1);ones(B,1)];
if d == 2
[betaActual(1,trial), linem(1,trial), lineb(1,trial)] = weakseparator(a,b,d);
else
[betaActual(1,trial), a1(1,trial),a2(1,trial),a3(1,trial),a4(1,trial)] = weakseparator(a,b,d);
end
betaActual(1,trial) = betaActual(1,trial)/(A+B);
U = [X Y];
%p0 = U(1,:);
%U = U([2:end],:);
C = dpptuples(U,d);
%C = dptuples(U,d,p0);
comblength = size(C,3);

separability = zeros(1,comblength);
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
    separability(1,i) = testsep(atest,btest,d);
    
end
separableTuples = sum(separability(1,:));


alpha = sum(separability)/comblength;
FK(1,trial) = alpha;
betaTheory(1,trial) = (rFinder(separableTuples,(A+B),d+1)+d+1+1)/(A+B);
if mod(trial,200) == 0
    save('RandFK.5_n_'+string(A+B)+'d_'+string(d)+'_T_'+string(trial));
end
if plotsingle == 1;
scatter(ahist(:,1,T),ahist(:,2,T)); hold on; scatter(bhist(:,1,T),bhist(:,2,T));hold on;
x = linspace(-20,20);
y = linem(T)*x+lineb(T);
plot(x,y);
%%%3dplotting
plot3(a(:,1),a(:,2),a(:,3),'o'); hold on; plot3(b(:,1),b(:,2),b(:,3),'x');
[x y] = meshgrid(-20:1:20);
hold on;
z = -a1(T)/a3(T).*x-a2(T)/a3(T).*y+a4(T)/a3(T);
surf(x,y,z);
end

%[r pk fk] = rFinder(alpha,(A+B),d);
end
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
C = zeros(4,size(U,2),Cnsize)
for i = 1:Cnsize
    for j = 1:d+1
    C(j,:,i) = U(Cn(i,j),:);
    %C(2,:,i) = U(Cn(i,2),:);
    %C(3,:,i) = U(Cn(i,3),:);
    %C(4,:,i) = p;
    end
    C(d+2,:,i) = p;
end
end
function [separable] = testsep(a,b,d)
sizea = size(a,1);
sizeb = size(b,1);
if sizea== d+2 || sizeb == d+2
    separable = 1;
   
else
A1 = [-a; b]
bs1 = [ones(1,sizea) -ones(1,sizeb)];
A2= [a;-b]
bs2 = [-ones(1,sizea) ones(1,sizeb)];
f = zeros(1,d);
separable1 = any(logical(linprog(f,A1,bs1)));
separable2 = any(logical(linprog(f,A2,bs2)));
separable = separable1||separable2;
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
function [a,b] = setBeta(A,B,beta)
n = A+B;
v = floor(beta*n/2);

for i = 1:v
a(i,:) = [3*cos(2*pi/v*i)-10 3*sin(2*pi/v*i)];
b(i,:) = [3*cos(2*pi/v*i)+10 3*sin(2*pi/v*i)];
end
for i = v+1:A
a(i,:) = [cos(2*pi/(A-v)*i)+10 sin(2*pi/(A-v)*i)];
b(i,:) = [cos(2*pi/(B-v)*i)-10 sin(2*pi/(B-v)*i)];
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
function rmax = rFinder(alphaN,n,d)
k = d+1;
%limit = alpha*nchoosek(n,k);
limit = alphaN;
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
