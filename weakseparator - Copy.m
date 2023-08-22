
function varargout = weakseparator(a,b,d)


U = [a;b];
maxdepth = 0;
n = size(U,1);
if d==2
    for i = 1:n-1
        for j = i+1:n
            x = [U(i,1) U(j,1)];
            y = [U(i,2) U(j,2)];
            p = polyfit(x,y,1);
            tvec = [-p(1) 1];
            temp = 0;
            for m = 1:n
                if (m ~= i && m~= j)
                    if m<= size(a,1)
                        if U(m,2)>=p(1,1)*U(m,1)+p(1,2)
                            temp= temp+1;
                        end
            
                    else
                        if U(m,2)<=p(1,1)*U(m,1)+p(1,2)
                            temp= temp+1;
                        end
                    end
                end
            end
    if temp+2>maxdepth
        maxdepth = temp+2;
        linem = p(1,1);
        lineb = p(1,2);
    end
    
    temp = 0;
    for m = 1:n
        if (m ~= i && m~= j)
            if m <= size(a,1);
                if U(m,2)<=p(1,1)*U(m,1)+p(1,2)
                    temp= temp+1;
                end
            else
                if U(m,2)>=p(1,1)*U(m,1)+p(1,2)
                    temp= temp+1;
                end 
            end
        end
    end
            if temp+2>maxdepth
                maxdepth = temp+2;
                linem = p(1);
                lineb = p(1,2);
            end
        end
    end
elseif d == 3
    Cn = combnk(1:n,3);
    Cnsize = size(Cn,1);
    for i = 1:Cnsize
       x = [U(Cn(i,1),1) U(Cn(i,2),1) U(Cn(i,3),1)];
       y = [U(Cn(i,1),2) U(Cn(i,2),2) U(Cn(i,3),2)];
       z = [U(Cn(i,1),3) U(Cn(i,2),3) U(Cn(i,3),3)];
       [a1,a2,a3,a4] = planefit([x;y;z]);
       temp1= 0; 
       temp2 =0;
    for m = 1:n
        if  ~any(Cn(i,:) == m)
            if m<= size(a,1)
                %Test Conditions;
                if a1*U(m,1)+a2*U(m,2)+a3*U(m,3) > a4
                    temp1 = temp1+1;
                else
                    temp2 = temp2+1;
                end
            else
               if a1*U(m,1)+a2*U(m,2)+a3*U(m,3)< a4
                   temp1 = temp1 +1;
               else
                   temp2 = temp2+1;
               end
            end
            if max(temp1,temp2)+3>maxdepth
                maxdepth = max(temp1,temp2)+3;
                aa1 = a1; aa2 = a2; aa3 = a3; aa4 = a4;
            end
        end
    end
    end
end
if nargout == 3
    varargout = {maxdepth, linem, lineb};
end
if nargout == 5
    varargout = {maxdepth, aa1, aa2, aa3, aa4};
end

function varargout=planefit(xyz)
%Fit 3D plane to given (x,y,z) data
%
%    [a,b,c,d] =  planefit(xyz)   
%     [abc,d]  =  planefit(xyz)
%      hcoeff  =  planefit(xyz)
%           
%IN:        
%           
%    xyz:   3xN input array with coordinates organized in columns
%           [x1,x2,x3...;y1,y2,y3,...;z1,z2,z3,...]
%                 
%OUT:       
%           
%    [a,b,c,d] : Coefficients of fitted plane equation a*x+b*y+c*z=d
%    [abc,d]   : Coefficients in 2-argument form where abc=[a,b,c]
%     hcoeff   : Homogeneous coefficient vector hcoeff=[a,b,c,-d]
 if size(xyz,1)~=3
    error 'Input xyz matrix must be 3xN' 
 end
 xyz=xyz.';
 
 mu=mean(xyz,1);
 
 [~,~,V]=svd(xyz-mu,0);
 
 normal=V(:,end).';
 cons=normal*mu';
 
 switch nargout
     
     case {0,1}
 
      varargout={[normal,-cons]};
      
     case 2
         
      varargout={ normal, cons};
      
     case 4
       
      varargout=[num2cell(normal),{cons}]; 
         
     otherwise
         
      error 'Must be 1,2, or 4 output args'
          
 end
 
end


end

