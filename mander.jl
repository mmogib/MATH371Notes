using Plots
using LinearAlgebra
Points = Vector{Tuple{Vector{Float64},Int}}(undef, 1000);

manderblut(c::Number) = z::Number -> z^2 + c
for i in 1:1000
    c1 = rand(-2:0.001:2)
    c2 = rand(-2:0.001:2)
    c = c1 + c2 * im
    g = manderblut(c)
    z = g(Complex(0))
    iters = 0
    while true
        iters += 1
        if iters > 100
            Points[i] = ([real(c), imag(c)], i)
            break
        end
        if norm(z) > 2
            Points[i] = ([real(c), imag(c)], i)
            break
        end
        z = g(z)

    end
end
x = map(x -> x[1][1], Points)
y = map(x -> x[1][2], Points)
clrs = map(x -> RGBA(0, 0, 0, x[2] / 1000), Points)
p = plot(xlims=(-10, 10), ylims=(-10, 10), label=nothing)
for i in 1:1000
    p = scatter(p, [x[i]], [y[i]], c=clrs[i], lw=0.2, label=nothing)
end
p = plot(p, label=nothing)