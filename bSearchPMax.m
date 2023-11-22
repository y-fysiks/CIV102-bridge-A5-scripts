lo = 1;
hi = 2000;
maxLoad = 0;
FOS = zeros(1, 8);

while lo < hi
    mid = ceil((lo + hi) / 2);
    [fails, maxLoad, FOS] = checkPFail(mid);
    if fails
        hi = mid - 1;
    else
        lo = mid;
    end
end

maxLoad
FOS