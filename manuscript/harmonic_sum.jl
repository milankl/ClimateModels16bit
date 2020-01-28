function harmonic_sum(::Type{T},steps::Int) where T

    s = zero(T)
    o = one(T)

    for i in 1:steps

        s_old = s
        s += o/T(i)

        if s == s_old    # check for convergence
            println(Float64(s),i)
            break
        end
    end
end
