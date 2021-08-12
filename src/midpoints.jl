

function midpoint(lo, hi, step)
    if isfinite(lo) && isfinite(hi)
        return (lo + hi)/2
    elseif isfinite(lo)
        return lo + step
    elseif isfinite(hi)
        return hi - step
    end
    0
end


shift(val, shift) = val + shift
