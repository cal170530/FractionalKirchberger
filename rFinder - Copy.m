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
