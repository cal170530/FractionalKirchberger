function varargout = FK_Hist(d,ahist,bhist,T,instances)
%FK_experiment - given two point classes of the same number of points in
%dimension d, this function: 
% 1) Produces linearly separable configurations of points(instances) in dimension d. 
% 2) For each instance produced in 1, a set of T new configurations are produced by swapping 1,2,..,T
%    pairs of labels between points of each class. This should produce a
%    range of beta values, where beta is the proportion of points correctly
%    classified by an optimal linear separator.
% ----The remaining steps are repeated for each subconfiguration(2) of each
% ---- point configuration(1):
%    3) Brute force search for the optimal linear classifier and obtain beta. 
%    4) Compute alpha. 
%    5) For a range of sample sizes over the (d+2) subsets, compute alpha_sample
%    6) from alpha and alpha_sample, obtain beta_theory and beta_sample. 
%    7) Store alpha, alpha_sample, beta, beta_theory, beta_sample.  Store
%    point configuration data. 
% 
%
%EXAMPLE:
% [betaActual,betaTheory,betaTheory_sample,betaTheory_sample2,alpha,alpha_sample,alpha_sample2] = FK_Experiment(2,10,5,1)
%  
switch nargin
    case(2)
        T = 10;
        instances = 10;
    case(3)
        instances = 10;
end

samplesteps = 10;
betaActual = zeros(instances,T);
betaTheory = zeros(instances,T);
betaTheory_sample = zeros(instances,samplesteps,T);
betaTheory_sample2 = zeros(instances,samplesteps,T);
alpha = zeros(instances,T);
alpha_sample = zeros(instances,T);
alpha_sample2 = zeros(instances,T);

cf = zeros(d+1,instances,T);



for trial = 1:T
    for inst = 1:instances
        a = sparse(ahist(:,:,:,trial));
        b = sparse(bhist(:,:,:,trial));


  

        A = size(a,1);
        B = size(b,1);

     



        X = [a;b];
        Y= [ zeros(A,1);ones(B,1)];
        U = [X Y];
        %[maxdepth, cf(:,inst,trial)] = weakseparator(a,b,d);

        %betaActual(inst,trial) = maxdepth/(A+B);

        for samp_percent = .005:.005:.05
            [C,Cn] = dptuples(U,d);
            comblength = size(C,3);
            separableTuples = 0;
            sampleSize= round((samp_percent)*comblength); 

%%%%%Setup for Sampling Local Test Results%%%%%%%%
            separability = zeros(1,comblength);
            sampleseparability = zeros(1,sampleSize);
            samples= datasample(1:comblength,sampleSize,'Replace',false);
            separableSamples = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Setup for Sampling point subsets%%%%%%%%%
            sampleSize2 = 0;
            for i = d+2:(A+B)
                if nchoosek(i,d+1)>sampleSize
                    sampleSize2 = i;
                    break;
                end
            end
            %sampleSize2 = round(samp_percent*(A+B));
            samplePoints2 = datasample(1:(A+B),sampleSize2,'Replace',false);
            sampleTuples2 = combnk(samplePoints2,d+1);
            separability2 = zeros(1,comblength);
            sampleseparability2 = zeros(1,size(sampleTuples2,1));
            Cn = sort(Cn,2);
            sampleTuples2 = sort(sampleTuples2,2);
            sampleTuples2 = sortrows(sampleTuples2);
            separableSamples2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for dir = -1:2:1
                    j=1;
                    k = 1;
    
                    for i = 1:sampleSize
                        aindex= find(C(:,d+1,samples(1,i)) == 0);
                        bindex = find(C(:,d+1,samples(1,i)) == 1);
                        if (aindex)
                            atest = C(aindex,1:d,i);
                        else
                            atest = [];
                        end
                        if(bindex)
                            btest = C(bindex,1:d,i);  
                        else
                            btest = [];
                        end
                        separability(1,i) = testsep(atest,btest,d,dir);
                        if any(samples(1,:) == i)
                            sampleseparability(1,j) = separability(1,i);
                            j = j+1;
                        end
                        %{
                        if any(ismember(Cn(i,:),sampleTuples2,'rows'))
                            sampleseparability2(1,k) = separability(1,i);
                            k = k+1;
                        end
                        %}
                    end
                    if (sum(separability(1,:))>separableTuples)
                        separableTuples = sum(separability(1,:));
                    end
                    if (sum(sampleseparability(1,:)) > separableSamples)
                        separableSamples = sum(sampleseparability(1,:));
                    end

                    %{
                    if (sum(sampleseparability2(1,:)) > separableSamples2)
                        separableSamples2 = sum(sampleseparability2(1,:));
                    end
                    %}
                end
            for dir = -1:2:1
                    j=1;
                    k = 1;
    
                    for i = 1:size(sampleTuples2,1)
                        aindex = [];
                        bindex = [];
                        for n = 1:d+1
                            if U(sampleTuples2(i,n),d+1)==0
                                aindex = [aindex; sampleTuples2(i,n)]
                            else
                                bindex = [bindex; sampleTuples2(i,n)]
                            end
                        end
                        if (aindex)
                            atest = U(aindex,1:d);
                        else
                            atest = [];
                        end
                        if(bindex)
                            btest = U(bindex,1:d);  
                        else
                            btest = [];
                        end
                        sampleseparability2(1,i) = testsep(atest,btest,d,dir);
                     
                     
                    end
                    
                    if (sum(sampleseparability2(1,:)) > separableSamples2)
                        separableSamples2 = sum(sampleseparability2(1,:));
                    end
                    
                end
            alpha(inst,trial) = separableTuples/comblength;
            alpha_sample(inst,trial) = separableSamples/sampleSize;
            alpha_sample2(inst,trial) = separableSamples2/size(sampleTuples2,1);
            betaTheory(inst,trial) = (rFinder(separableTuples,(A+B),d)+d+1)/(A+B);
            betaTheory_sample(inst,trial) = (rFinder(alpha_sample(1,trial)*comblength,(A+B),d)+d+1)/(A+B);
            betaTheory_sample2(inst,trial) = (rFinder(alpha_sample2(1,trial)*comblength,(A+B),d)+d+1)/(A+B);
            if trial == 10
                save('FK_Sampling_Global_Local_T_'+string(trial));
            end
        end 

    end
end
switch nargout
    case(3)
        varargout = {betaActual,betaTheory,betaTheory_sample};
    case(4)
        varargout = {betaActual,betaTheory,betaTheory_sample,alpha};
    case(5)
        varargout = {betaActual,betaTheory,betaTheory_sample,alpha,alpha_sample};
    case(7)
        varargout= {betaActual,betaTheory,betaTheory_sample,betaTheory_sample2,alpha,alpha_sample,alpha_sample2}
   
end



function [C] = dpptuples(U,d)
length = size(U,1);   
Cn = combnk(1:length,d+2);
Cnsize = size(Cn,1);
C = zeros(4,size(U,2),Cnsize)
    for i = 1:Cnsize
        for j = 1:d+2
        C(j,:,i) = U(Cn(i,j),:);
        end
    end
end

function[C,Cn] = dptuples(U,d,p)
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
        A1 = [-a; b];
        bs1 = [ones(1,sizea) -ones(1,sizeb)];
        A2= [a;-b];
        bs2 = [-ones(1,sizea) ones(1,sizeb)];
        f = zeros(1,d);
        if dir >0
            separable = any(logical(linprog(f,A1,bs1)));
        else
            separable = any(logical(linprog(f,A2,bs2)));
        end

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
function [a,b] = dividePoints(n,d)
    a = zeros(n,d);
    b = zeros(n,d);
    c = zeros(1,d);
    for i = 1:d-1
        c(1,i) = randi(5)/randi(5);
    end
    c(1,d) = randi(10)-6;

    for i = 1:n
        for j = 1:d-1
            a(i,j) = (rand()-.5)*50;
            b(i,j) = (rand()-.5)*50;
            a(i,d)= a(i,d)+c(1,j)*a(i,j);
            b(i,d) = b(i,d)+c(1,j)*b(i,j);
        end
        a(i,d) = a(i,d)+c(1,d);
        b(i,d) = b(i,d)+c(1,d);
        a(i,d) = a(i,d)+(50*(d-1)*5-a(i,d))*rand();
        b(i,d) = b(i,d)-(50*(d-1)*5-b(i,d))*rand();
    end

end

end