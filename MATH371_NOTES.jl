### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 65bdc140-2f92-11ef-1cbe-31065d820068
begin
    using CommonMark
    using PlutoUI, PlutoExtras
    using Plots, PlotThemes, LaTeXStrings
    using Latexify
    using HypertextLiteral
    using Colors
    using LinearAlgebra, Random, Printf
    using Symbolics
    using QRCoders
    using PrettyTables
    using NonlinearSolve
    # using ForwardDiff
    using Integrals
	using OrdinaryDiffEq
end

# ╔═╡ 8ca0d1c5-166d-44f0-a17e-a6207c19459a
initialize_eqref()

# ╔═╡ dbc554cb-3061-48b8-aaa2-c6c13e7a3d75
@htl("""
<style>
@import url("https://mmogib.github.io/math102/custom.css");


</style>
""")

# ╔═╡ 21331300-b8cb-410f-88a3-1feb63e55659
TableOfContents(title="MATH371")

# ╔═╡ a279f75a-1d5a-4c0e-aae6-a168dda0a277
begin

    exportqrcode("https://mshahrani.website/teaching/", "website.png"; width=1)

end

# ╔═╡ f581c560-afc6-460d-bb7e-48e13c1909d4
begin
    struct LocalImage
        filename
    end

    function Base.show(io::IO, ::MIME"image/png", w::LocalImage)
        write(io, read(w.filename))
    end
end

# ╔═╡ 7d454fff-b638-4678-92f6-12d816b541b1
LocalImage("./website.png")

# ╔═╡ 830e0fb0-8bfa-48b6-9210-cc0dc009d042
md"__My Website__"

# ╔═╡ d7267eec-8b39-454e-b1ec-c7ae580190c4
md"# 1.1: Review of Calculus: Taylor Polynomials and Series"

# ╔═╡ cb26e993-d7c1-4e69-82f3-dcb20d1a4f37
begin
    @syms x::Real
    Df = Differential(x)
    D(f, x) = begin
        expand_derivatives(Df(f(x)))
    end
    df(f, n, x0) = begin
        val = if n == 0
            f(x)
        else
            reduce((acc, val) -> begin
                    D(t -> substitute(acc, Dict(x => t)), x)
                end, 2:n; init=D(f, x))
        end
        substitute(val, Dict(x => x0))
    end
end

# ╔═╡ b54035ab-2813-4cdc-887a-e11625ede4aa
slider1h = @bind slider1 Slider(0:20, show_value=true);
"";

# ╔═╡ 9e344be4-e96a-4b1f-9bac-52deb66c2658


# ╔═╡ d60047cf-47fd-406f-b100-be2b1d195612
md"n = $slider1h"

# ╔═╡ 7e6509af-4072-45f4-a3cf-d2891c586d58
L"$P_{%$slider1}(x)=$"

# ╔═╡ df80d1bc-1285-48f3-8de1-5bf9631f078c
L"$R_{%$slider1}(x)=$"

# ╔═╡ 2d460985-77a6-4453-b639-3312c98b742f
md"When ``x_0=0.01``"

# ╔═╡ a4214d60-b66e-4819-84ad-37653e45184e
cm"""
Previous Example illustrates __the two objectives of numerical analysis__:
- (i) Find an approximation to the solution of a given problem.
- (ii) Determine a bound for the accuracy of the approximation.
"""

# ╔═╡ f610361d-e308-49be-8100-ad8b3882bac0
md"# 1.2 Round-off Errors and Computer Arithmetic"

# ╔═╡ e4a92b5e-18a7-4007-aab3-3eeca447b8de
md"""
## Binary Machine Numbers
- In 1985, the IEEE (Institute for Electrical and Electronic Engineers) published a report called Binary Floating Point Arithmetic Standard [754-1985](https://www.dropbox.com/scl/fi/xa1ikgu16lcv3bo9f2oz9/754-1985.pdf?rlkey=l1vn0cbev2ss23pu3m443t5je&raw= 1). 
- This was updated in 2008 as IEEE 754-2008.
- The current version, [IEEE 754-2019](https://www.dropbox.com/scl/fi/vggm34gpdlvr2x3s87gy8/754-2019.pdf?rlkey=xbjgqy23exr6q91jj2rwzrut9&raw=1) published in July 2019, is derived from and replaces IEEE 754-2008, [READ MORE](https://754r.ucbtest.org/)

> This provides standards for 
> - __arithmetic formats__: sets of binary and decimal floating-point data, which consist of finite numbers (including signed zeros and subnormal numbers), infinities, and special "not a number" values (NaNs)
> - __interchange formats__: encodings (bit strings) that may be used to exchange floating-point data in an efficient and compact form
> - __rounding rules__: properties to be satisfied when rounding numbers during arithmetic and conversions
> - __operations__: arithmetic and other operations (such as trigonometric functions) on arithmetic formats
> - __exception handling__: indications of exceptional conditions (such as division by zero, overflow, etc.)

"""

# ╔═╡ 11c16f32-f86b-43ba-8a58-532609c226e0
cm"""
__A 64-bit (binary digit) representation is used for a real number.__
- The first bit is a __sign indicator__, denoted ``s``. 
- This is followed by an 11-bit exponent, ``c``, called the __characteristic__, and - a 52 -bit binary fraction, ``f``, called the __mantissa__. 

The base for the exponent is 2 .

> To save storage and provide a unique representation for each floating-point number, a normalization is imposed. Using this system gives a floating-point number of the form
```math
(-1)^s 2^{c-1023}(1+f)
```

"""


# ╔═╡ c9819b7f-4bf8-46e3-a4bc-e991bf26ce7a


# ╔═╡ ed835d55-deaf-44b5-b360-f8eff7800f6d
cm"""
- What is the next smallest machine real number to the previous number ?
- What is the next largest machine real number to the previous number ?
"""

# ╔═╡ 36211411-8f83-4373-8dd0-7b7f334f8bc1
cm"""
- What is the smallest real number can be represented ?
- What is the largest real number can be represented?

"""

# ╔═╡ 66fa9c4b-f2f3-44a1-8755-bc13ffccaba1
let
    x = split("27.56640625", '.')
    # sum(float(Int(x[2][i])//(2^(i))) for i in 1:length(x[2]))
    bs = bitstring(parse(Float64, "27.56640625"))
    reinterpret(Float64, parse(Int, "0100000000111011100100001111111111111111111111111111111111111111", base=2))

end

# ╔═╡ ff2fe07a-d35b-4d5b-bf95-f669fafc2132
md"## Decimal Machine Numbers"

# ╔═╡ 9ddc2701-bf18-4f96-9394-36a95b399252
cm"""
We assume that machine numbers are represented in the normalized decimal floating-point form
```math
\pm 0 . d_1 d_2 \ldots d_k \times 10^n, \quad 1 \leq d_1 \leq 9, \quad \text { and } \quad 0 \leq d_i \leq 9,
```
for each ``i=2, \ldots, k``. Numbers of this form are called ``k``-digit decimal machine numbers.
- Any positive real number within the numerical range of the machine can be normalized to the form
```math
y=0 . d_1 d_2 \ldots d_k d_{k+1} d_{k+2} \ldots \times 10^n
```
- BUT our computer has to terminate the __mantissa__ at ``k`` decimal digits.
- The floating-point form of ``y``, denoted ``\text{fl}(y)``, is obtained by terminating the mantissa of ``y`` at ``k`` decimal digits. 
- There are two common ways of performing this termination. 
  - __One method__, called __chopping__, is to simply chop off the digits ``d_{k+1} d_{k+2} \ldots``. This produces the floating-point form
```math
\text{fl}(y)=0 . d_1 d_2 \ldots d_k \times 10^n .
```
-
  - The other method, called __rounding__, adds ``5 \times 10^{n-(k+1)}`` to ``y`` and then chops the result to obtain a number of the form
```math
\text{fl}(y)=0 . \delta_1 \delta_2 \ldots \delta_k \times 10^n .
```

> For rounding, when ``d_{k+1} \geq 5``, we add 1 to ``d_k`` to obtain ``\text{fl}(y)``; that is, we round up. When ``d_{k+1}<5``, we simply chop off all but the first ``k`` digits; that is, round down. If we round down, then ``\delta_i=d_i``, for each ``i=1,2, \ldots, k``. However, if we round up, the digits (and even the exponent) might change.
"""

# ╔═╡ 7caaa18e-8ef8-4ae5-bd1c-e9fd0536b901
π

# ╔═╡ 43e1134e-020c-4fff-8164-050f0fcb9c38
# let 
# 	n = 0.31415926535897e1
# 	# chooping
# 	n_chopped =0.31415e1 
# 	n_rounded = 0.31416e1
# end

# ╔═╡ 8c6005c5-e947-4a8a-8083-78ab5f84f07b
# let
# 	function actualErro(p,pstar)
# 		p-pstar
# 	end
# 	function absoluteErro(p,pstar)
# 		abs(p-pstar)
# 	end
# 	function relativeErro(p,pstar)
# 		abs(p-pstar)/abs(p)
# 	end
# 	p,pstar=0.3000e1,0.3100e1
# 	actualErro(p,pstar)
# 	absoluteErro(p,pstar)
# 	relativeErro(p,pstar)
# end

# ╔═╡ c8d040f0-3264-42c6-84c5-be8c633c36e4
let
    # (a)
    aerror1 = 0.1e2
end


# ╔═╡ 786351e5-47cc-4d61-a47d-771bc47aa3a9
cm"""
For machine decimal representations of numbers we have
```math
\left|\frac{y-f l(y)}{y}\right| \text {. }
```

If ``k`` decimal digits and chopping are used for the machine representation of
```math
y=0 . d_1 d_2 \ldots d_k d_{k+1} \ldots \times 10^n
```
then
```math
\begin{aligned}
\left|\frac{y-f l(y)}{y}\right| & =\left|\frac{0 . d_1 d_2 \ldots d_k d_{k+1} \ldots \times 10^n-0 . d_1 d_2 \ldots d_k \times 10^n}{0 . d_1 d_2 \ldots \times 10^n}\right| \\
& =\left|\frac{0 . d_{k+1} d_{k+2} \ldots \times 10^{n-k}}{0 . d_1 d_2 \ldots \times 10^n}\right|=\left|\frac{0 . d_{k+1} d_{k+2} \ldots}{0 . d_1 d_2 \ldots}\right| \times 10^{-k} .
\end{aligned}
```

Since ``d_1 \neq 0``, the minimal value of the denominator is 0.1 . The numerator is bounded above by 1 . As a consequence,
```math
\left|\frac{y-f l(y)}{y}\right| \leq \frac{1}{0.1} \times 10^{-k}=10^{-k+1} .
```

"""

# ╔═╡ 2e3eefec-a29b-4dff-a863-eabbb4a3c7ef
md"## Finite-Digit Arithmetic"

# ╔═╡ 712e7cbd-1280-4885-8adb-d83de2c26560
begin
    struct FloatK
        m::Int
        n::Int
        s::Int
        digits::Int
        rounding::String
    end
    Base.convert(::Type{Float64}, x::FloatK) = x.s * x.m * (10.0)^(x.n - x.digits)
    Base.show(io::IO, ::MIME"text/plain", n::FloatK) = print(io, n.s == 1 ? "" : "-", " 0.", n.m, "× 10^", n.n)
    Base.show(io::IO, ::MIME"application/x-tex", n::FloatK) = print(io, n.s == 1 ? "" : "-", "0.", n.m, "× 10^", n.n)
    FloatK(n::Float64, dgts::Int, rounding::String) = begin
        s = n >= 0 ? 1 : -1
        fpart, ipart = modf(n)
        ipart = abs(ipart)
        fpart = abs(fpart)
        no_ipart = (ipart == 0.0)
        fparts = "$fpart"
        fpartstr, expntn = if 'e' in fparts
            fe = findfirst('e', fparts)
            np = fparts[1:1] * fparts[3:fe-1]
            nex = parse(Int, fparts[fe+1:end]) + 1
            fraction_str, fraction_exp = if no_ipart
                repeat('0', abs(nex)) * np, nex
            else
                np, length("$(Int(ipart))")
            end
            fraction_str, fraction_exp
        else
            np = fparts[3:end]
            nex = findfirst(x -> x != '0', np)
            fraction_str, fraction_exp = if no_ipart
                np[nex:end], -nex + 1
            else
                np, length("$(Int(ipart))")
            end
            fraction_str, fraction_exp
        end
        fpartstr = fpartstr * repeat('0', 5 * dgts)
        iparttstr = parse(Int, split("$ipart", ".")[1])
        glued = no_ipart ? fpartstr : "$iparttstr" * fpartstr


        nstr = if rounding == "chop"
            glued[1:dgts]
        else
            glued_rounded = replace("$(round(parse(Int, glued[2:dgts+3]) * 10.0^(-dgts - 2), digits=dgts))", '.' => "") * repeat('0', 2 * dgts)
            last_digit = parse(Int, glued_rounded[dgts+1])
            correction = last_digit >= 5 ? 1 : 0
            temp = "$(parse(Int, glued[1:dgts]) + correction))"
            temp[1:dgts]
        end

        m = parse(Int, nstr)
        FloatK(m, expntn, s, dgts, rounding)
    end
    FloatK(x::FloatK) = FloatK(x.m, x.n, x.s, x.digits, x.rounding)
    FloatK(x::FloatK, rounding::String) = FloatK(x.m, x.n, x.s, x.digits, rounding)
    FloatK(x::FloatK, digits::Int, rounding::String) = FloatK(x.m, x.n, x.s, digits, rounding)
    FloatK(n::Float64, dgts::Int) = FloatK(n, dgts, "chop")
    FloatK(n::Float64, rounding::String) = FloatK(n, 5, rounding)
    FloatK(n::Float64) = FloatK(n, 5, "chop")

    Base.:+(a::FloatK, b::FloatK) = FloatK(convert(Float64, a) + convert(Float64, b), a.digits, a.rounding)
    Base.:+(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, a) + b, a.digits, a.rounding)
    Base.:+(b::T where {T<:Real}, a::FloatK) = FloatK(a, b)

    Base.:-(a::FloatK) = FloatK(-convert(Float64, a))
    Base.:-(a::FloatK, b::FloatK) = FloatK(convert(Float64, a) - convert(Float64, b), a.digits, a.rounding)
    Base.:-(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, a) - b, a.digits, a.rounding)
    Base.:-(b::T where {T<:Real}, a::FloatK) = FloatK(b - convert(Float64, a), a.digits, a.rounding)

    Base.:*(a::FloatK, b::FloatK) = FloatK(convert(Float64, a) * convert(Float64, b), a.digits, a.rounding)
    Base.:*(a::T where {T<:Real}, b::FloatK) = FloatK(a * convert(Float64, b), b.digits, b.rounding)
    Base.:*(a::FloatK, b::T where {T<:Real}) = FloatK(b, a)

    Base.:^(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, a)^b, a.digits, a.rounding)
    Base.:sqrt(a::FloatK) = FloatK(sqrt(convert(Float64, a)), a.digits, a.rounding)
    Base.:÷(a::FloatK, b::FloatK) = FloatK(convert(Float64, a) / convert(Float64, b), a.digits, a.rounding)
    Base.:÷(a::T where {T<:Real}, b::FloatK) = FloatK(a / convert(Float64, b), b.digits, b.rounding)
    Base.:÷(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, b) / a, a.digits, a.rounding)

    Base.:/(a::FloatK, b::FloatK) = FloatK(convert(Float64, a) / convert(Float64, b), a.digits, a.rounding)
    Base.:/(a::T where {T<:Real}, b::FloatK) = FloatK(a / convert(Float64, b), b.digits, b.rounding)
    Base.:/(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, b) / a, a.digits, a.rounding)

    function createFiniteDigitSystem(; digits::Int=5, truncation::String="chop")
        fl(x) = FloatK(x, digits, truncation)
        ⊕(a::Float64, b::Float64) = FloatK(FloatK(a, digits, truncation) + FloatK(b, digits, truncation), digits, truncation)
        ⊖(a::Float64, b::Float64) = FloatK(FloatK(a, digits, truncation) - FloatK(b, digits, truncation), digits, truncation)
        ⊗(a::Float64, b::Float64) = FloatK(FloatK(a, digits, truncation) * FloatK(b, digits, truncation), digits, truncation)
        ⨸(a::Float64, b::Float64) = FloatK(FloatK(a, digits, truncation) ÷ FloatK(b, digits, truncation), digits, truncation)
        fl, ⊕, ⊖, ⊗, ⨸
    end

    function re(p::T where {T<:Number}, ps::FloatK)
        abs(convert(Float64, (p - ps)) / p)
    end
end

# ╔═╡ e223e89c-af0c-44ed-bb83-fd51124c7899
let
    f(x) = cos(x)
    P2(x) = 1 - 0.5 * x^2
    x = 3
    f(x), P2(x)
end

# ╔═╡ ee13d4c2-edce-4d7c-8ace-103aab1f7ba0
let
    f(x) = cos(x)
    x0 = 0
    # slider1

    P(n) = x0 -> sum(((x - x0)^i) * df(f, i, x0) / factorial(i) for i in 0:n)
    Pn = P(slider1)(x0)
    Pn

end

# ╔═╡ afbd4725-9c76-4453-a62a-ebac185655a2
let
    @syms ζ::Real
    x0 = 0
    n = slider1
    f(x) = cos(x)
    ((x - x0)^(n + 1)) * df(f, n + 1, ζ) / factorial(n + 1)
end

# ╔═╡ 422a3ccc-1e62-407b-bd3e-0703f1a5a33d
let
    f(x) = cos(x)
    x0 = 0
    P(n) = x0 -> sum(((x - x0)^i) * df(f, i, x0) / factorial(i) for i in 0:n)
    Pn = P(slider1)(x0)
    p1 = plot(f; framestyle=:origin, label=L"f(x)", line=(1, :red))
    labels = [L"P_{%$i}(x)" for i in 0:2:slider1]
    plot(p1, [t -> substitute(P(i)(x0), Dict(x => t)) for i in 0:2:slider1];
        label=reshape(labels, 1, length(labels)),
        xlimits=(-12, 12), ylimits=(-1, 1.5)
    )
end

# ╔═╡ 9b0cd30c-3162-4423-b32f-88c9a01c1cbd
let
    x1 = 0.01
    n = slider1
    f(x) = cos(x)
    x0 = 0
    P(n) = x0 -> sum(((x - x0)^i) * df(f, i, x0) / factorial(i) for i in 0:n)
    vals = substitute(P(n)(x0), Dict(x => x1)), f(x1)
    Rn = ((x - x0)^(n + 1)) / factorial(n + 1) |> y -> substitute(y, Dict(x => x1))
    pntex = texeq("""
    P_{$n}($x1) = $(vals[1])
    	""")
    ftex = texeq("f($x1) = $(vals[2])")
    error = texeq("|\\textrm{Error}| \\leq $(Rn)")
    md"""
    $ftex

    $pntex

    $error
    """

end

# ╔═╡ 98652338-3d5c-40f9-9220-c4f6e94f532c
sqrt(2)

# ╔═╡ 2563dd79-3e17-4b81-82f3-43d6aea85b71
begin
    function bint64_to_dec(str::String)
        if length(str) > 64
            error("The binary string is too long. It should be less than 64.")
        end
        if length(filter(x -> x in ['0', '1'], str)) != length(str)
            error("The string is not valid.")
        end
        reduce((c, v) -> v[2] == '0' ? c : c + 2^(Float64(v[1] - 1)), zip(length(str):-1:1, str); init=0)
    end
    function float64_to_dec(str::String)
        if length(str) > 64
            error("The binary string is too long. It should be less than 64.")
        end
        if length(filter(x -> x in ['0', '1'], str)) != length(str)
            error("The string is not valid.")
        end
        lstr = length(str)
        str = lstr < 63 ? repeat('0', 64 - length(str)) * str : str
        # m,c,s= str[13:64], str[2:12], str[1]
        # mv = reduce((c,v)->v[2]=='0' ? c : c+1/2^(BigFloat(v[1])),zip(1:length(m),m);init=0)
        # cv = bint64_to_dec(c)
        # sv = s=='0' ? 1 : -1
        # # reduce((c,v)->v[2]=='0' ? c : c+2^(v[1]-1),zip(length(str):-1:1,str);init=0)

        # sv*2^(float(cv-1023))*(1+mv)
        reinterpret(Float64, parse(Int, str, base=2))
    end
end

# ╔═╡ 726a42a1-9746-4083-88f3-0661cd7396b5
bint64_to_dec("11111111111") - 1023

# ╔═╡ 4d8df2cc-74b1-445e-9cab-cb58af4d56bb
begin
    strs = "0100000000111011100100010000000000000000000000000000000000000000"
    float64_to_dec(strs)
    # 1+sum([1/2 1/8 1/16 1/32 1/256 1/4096])
end


# ╔═╡ 13ac7d12-0740-45df-831c-30670afd9f40
let
    a = Int(2)
    bint64_to_dec("11111111111")
    # parse(Int,bitstring(a);base=2)
    # typeof(str)
end

# ╔═╡ 239e7eaa-b8da-47d4-b642-8415ef1e0683
let
    binstr = "0100000000111011100100010000000000000000000000000000000000000000"
    next_small = "0100000000111011100100001111111111111111111111111111111111111111"
    next_largest = "0100000000111011100100010000000000000000000000000000000000000001"
    pn = float64_to_dec(next_small)
    n = float64_to_dec(binstr)
    nn = float64_to_dec.(next_largest)
    # 0.5(nn-n)
end

# ╔═╡ bb4ed82e-3b7e-4898-987f-9a9bac49c205
let
    smallest_number = String(repeat('0', 64))
    # smll_normalized =String(repeat('0',11)*repeat('0',52))
    # bitstring(1-1023)
    largest_number = '0' * repeat('1', 62) * '0'
    float64_to_dec(largest_number)
end

# ╔═╡ 75bac38d-3c93-4cad-8302-6eb60e40038b
cm"""
__BE CAREFULL__

```math
\large
0.1 + 0.2
```
$(0.1+0.2)
"""

# ╔═╡ d38e3f95-71c4-45dc-8ab9-0a216de91311
let
    strs = join(map(x -> "$(x*5*10^(-4.0))", [0.1 0.5 100 1000 5000 9990 10000]), ",")
    strs = L"%$strs"

end

# ╔═╡ 4db7b813-d0da-4fc6-8447-aae94fbae670
cm"""
In addition to inaccurate representation of numbers, the arithmetic performed in a computer is not exact. The arithmetic involves manipulating binary digits by various shifting, or logical, operations. Since the actual mechanics of these operations are not pertinent to this presentation, we shall devise our own approximation to computer arithmetic. Although our arithmetic will not give the exact picture, it suffices to explain the problems that occur.

We will assume a finite-digit arithmetic given by
```math
\begin{array}{ll}
x \oplus y=f l(f l(x)+f l(y)), & x \otimes y=f l(f l(x) \times f l(y)), \\
x \ominus y=f l(f l(x)-f l(y)), & x \circledast y=f l(f l(x) \div f l(y)) .
\end{array}
```
"""

# ╔═╡ 9010d6bd-26d2-4211-ab12-ac58d77edbbc
let
    re(p::T, ps::S) where {T<:Number,S<:Number} = abs((p - ps) / p)
    fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=5, truncation="chop")
    x = 5 / 7
    y = 1 / 3
    p, ps = x + y, x ⊕ y
    re(p, convert(Float64, ps))
    # x ⊖ y
    # x ⊗ y
    # x ⨸ y
end

# ╔═╡ b50497e4-e0f9-4da6-bba3-48628cee0173
let
    re(p::T, ps::S) where {T<:Number,S<:Number} = abs((p - ps) / p)
    fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=5, truncation="chop")
    x = 5 / 7
    y = 1 / 3
    u = 0.714251
    v = 98765.9
    w = 0.111111e-4
    p1 = 0.0003
    p2 = p1 ⊗ v
end

# ╔═╡ 379608bf-906d-45db-9cab-018a50fcbd07
let
    f(x) = x^2 + 62.10x + 1
    x1(a, b, c) = ((-b + sqrt(b^2 - 4 * a * c)) / 2a)
    x2(a, b, c) = ((-b - sqrt(b^2 - 4 * a * c)) / 2a)
    a, b, c = 1, 62.10, 1
    x1(a, b, c)
    x2(a, b, c)
    aa = cc = 1
    ba = 62.10
    ab2 = 0.3856e4
    # ab2= ba^2
    a4ac = 4
    asqr = 0.6206e2
    # asqr=sqrt(ab2-a4ac)
    anum1 = -0.0400
    ax1 = anum1 / 2.0
    anum2 = -0.1242e3
    ax2 = anum2 / 2.0


end

# ╔═╡ 1922d625-ba5e-4ea2-b050-a3c64077c742
let
    fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(digits=5, truncation="round")
    f(x) = x^2 + 62.10x + 1
    x1(a, b, c) = ((-b + sqrt(b^2 - 4 * a * c)) / 2a)
    x1f(a, b, c) = ((-fl(b) + fl(sqrt(fl(fl(fl(b)^2) - fl(4 * fl(a) * fl(c)))))) / fl(2a))
    x2(a, b, c) = ((-b - sqrt(b^2 - 4 * a * c)) / 2a)
    x2f(a, b, c) = ((-fl(b) - fl(sqrt(fl(fl(fl(b)^2) - fl(4 * fl(a) * fl(c)))))) / fl(2a))
    a, b, c = 1.0, 62.10, 1.0
    x1(a, b, c)
    x2f(a, b, c)
    # x2(a,b,c)
    # x2(fl(a),fl(b),fl(c))
end

# ╔═╡ 5ed7d51c-b5fa-4953-b2b4-c2040edced33
md"## Nested Arithmetic"

# ╔═╡ 987b5101-f2d6-4d03-98de-595bf6710d96
let
    fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=3, truncation="chop")
    f(x) = x^3 - 6.1 * x^2 + 3.2 * x + 1.5
    fn(x) = x * (x * (x - 6.1) + 3.2) + 1.5
    x = 4.71
    y = f(x)
    ys = f(fl(x))
    z = fn(x)
    zs = fn(fl(x))
    # fl(x^2,3,"round")
end

# ╔═╡ a7cc3418-4607-4e03-ab87-bab26530cc53
md"# 1.3 Algorithms and Convergence"

# ╔═╡ 78fc998d-a2c8-448b-90f5-d822e3513e8b
md"## Characterizing Algorithms"

# ╔═╡ 8cfd5ff2-cbf6-469e-af5a-050ecf2f3198
cm"""
| ``n`` | Computed ``\hat{p}_n`` | Correct ``p_n`` | Relative error |
| :--- | ---: | :--- | :--- |
| 0 | ``0.10000 \times 10^1`` | ``0.10000 \times 10^1`` |  |
| 1 | ``0.33333 \times 10^0`` | ``0.33333 \times 10^0`` |  |
| 2 | ``0.11110 \times 10^0`` | ``0.11111 \times 10^0`` | ``9 \times 10^{-5}`` |
| 3 | ``0.37000 \times 10^{-1}`` | ``0.37037 \times 10^{-1}`` | ``1 \times 10^{-3}`` |
| 4 | ``0.12230 \times 10^{-1}`` | ``0.12346 \times 10^{-1}`` | ``9 \times 10^{-3}`` |
| 5 | ``0.37660 \times 10^{-2}`` | ``0.41152 \times 10^{-2}`` | ``8 \times 10^{-2}`` |
| 6 | ``0.32300 \times 10^{-3}`` | ``0.13717 \times 10^{-2}`` | ``8 \times 10^{-1}`` |
| 7 | ``-0.26893 \times 10^{-2}`` | ``0.45725 \times 10^{-3}`` | ``7 \times 10^0`` |
| 8 | ``-0.92872 \times 10^{-2}`` | ``0.15242 \times 10^{-3}`` | ``6 \times 10^1`` |

__The error grows exponentially__
"""

# ╔═╡ 95e1b0f1-fec8-4c35-a527-a5351ed38d92
cm"""
| ``n`` | Computed ``\hat{p}_n`` | Correct ``p_n`` | Relative error |
| :--- | ---: | ---: | :---: |
| 0 | ``0.10000 \times 10^1`` | ``0.10000 \times 10^1`` |  |
| 1 | ``0.33333 \times 10^0`` | ``0.33333 \times 10^0`` |  |
| 2 | ``-0.33330 \times 10^0`` | ``-0.33333 \times 10^0`` | ``9 \times 10^{-5}`` |
| 3 | ``-0.10000 \times 10^1`` | ``-0.10000 \times 10^1`` | 0 |
| 4 | ``-0.16667 \times 10^1`` | ``-0.16667 \times 10^1`` | 0 |
| 5 | ``-0.23334 \times 10^1`` | ``-0.23333 \times 10^1`` | ``4 \times 10^{-5}`` |
| 6 | ``-0.30000 \times 10^1`` | ``-0.30000 \times 10^1`` | 0 |
| 7 | ``-0.36667 \times 10^1`` | ``-0.36667 \times 10^1`` | 0 |
| 8 | ``-0.43334 \times 10^1`` | ``-0.43333 \times 10^1`` | ``2 \times 10^{-5}`` |

__Error grows linearly__
"""

# ╔═╡ fdae435c-4d85-453e-a9ea-8383d31c5fe5
md"## Rates of Convergence"

# ╔═╡ f50aff72-d978-4c2a-8683-8127a14a4ea9
let
    @syms h::Real ζ::Real
    P(n) = sum(iseven(i) ? (-1)^i * h^(i) / (factorial(i)) : 0 for i in 0:n) + (iseven(n) ? sin(ζ) * h^(n + 1) / factorial(n + 1) : cos(ζ) * h^(n + 1) / factorial(n + 1))
    p3 = P(3)

end

# ╔═╡ 4dd54f9d-eb3f-49b4-9d56-41c5397ba001
md"# 2.1 The Bisection Method"

# ╔═╡ ee859bf9-4d20-46f5-9616-c932962cdbe2
function solveit(g, u0)
    f(u, p) = g(u)
    problem = NonlinearProblem(f, u0)
    sol = solve(problem)
    if sol.retcode == ReturnCode.Success
        sol.u, sol.stats, sol.alg
    else
        nothing
    end
end

# ╔═╡ 7d947cce-355f-4168-a76d-df5929d03be5
cm"""
__Bisection__

To find a solution to `` f(x) = 0 `` given the continuous function `` f `` on the interval ``[a, b]``, where `` f(a) `` and `` f(b) `` have opposite signs:

``\textbf{INPUT}`` endpoints `` a ``, `` b ``; tolerance ``\texttt{TOL}``; maximum number of iterations `` N_0 ``.

``\textbf{OUTPUT}`` approximate solution `` p `` or message of failure.

- ``\textbf{Step 1}`` Set `` i = 1 ``; `` FA = f(a). ``
- ``\textbf{Step 2}`` While `` i \leq N_0 `` do Steps 3--6.
- ``\textbf{Step 3}`` Set `` p = a + \frac{b - a}{2}; \quad \text{(Compute } p_i.\text{)} ``
        `` FP = f(p). ``
- ``\textbf{Step 4}`` If `` FP = 0 `` or `` \frac{b - a}{2} < \texttt{TOL} `` then 
        `` \text{OUTPUT } (p); \quad \text{(Procedure completed successfully.)} ``
        ``\text{STOP.}``

- ``\textbf{Step 5}`` Set `` i = i + 1 ``.
- ``\textbf{Step 6}`` `` \text{If }  FA \cdot FP > 0 `` then set `` a = p ``;  (Compute  ``a_i, b_i``.) `` FA = FP `` else set  `` b = p ``.  (FA is unchanged.)
- ``\textbf{Step 7}`` ``\text{OUTPUT}`` Method failed after ``N_0`` iterations, `` N_0 =``, ``N_0``);
    `` \text{(The procedure was unsuccessful.)} ``
    ``\text{STOP.}``
		"""

# ╔═╡ 9e6ed715-6bf8-4e60-b77e-f4f8e2118f02
begin
    function bisect(f, a, b, TOL, N0)
        i = 1
        FA = f(a)
        T = Matrix{Number}(undef, N0, 5)
        while i <= N0
            p = a + (b - a) / 2
            FP = f(p)
            T[i, :] = vcat(i, a, b, p, FP)
            if FP == 0 || ((b - a) / 2) < TOL
                TT = T[1:i, :]
                return p, TT
            end
            i = i + 1
            if FA * FP > 0
                a = p
                FA = FP
            else
                b = p
            end
        end
        @error("Maximum number of iterations reached")
    end
end

# ╔═╡ 3b64007d-762c-4bc5-8751-81ccb69ef376
let
    f(x) = x^2 * (x + 4) - 10
    # f(1),f(2)
    result = bisect(f, 1, 2, 1e-12, 50)
    p, T = result
    # pretty_table
    pretty_table(HTML, T[1:13, :], header=["n", "an", "bn", "pn", "f(pn)"])
    # f(p)
    # p
    # u, = solveit(f,2)
    # u

end

# ╔═╡ 505497ed-7060-48fe-ba2d-69a31413c267
md"# 2.2 Fixed-Point Iteration"

# ╔═╡ 1b28028a-92df-4d7a-942f-11ab9f4d06a7
begin
    function fixed_point(g, x0, ϵ; maxiters=50)
        n = maxiters
        # r(n,x0)=reduce((c,_)->abs(c-g(c)) <ϵ ? c : g(c),1:n,init=x0)
        xs = Vector{Float64}(undef, n)
        xs[1] = x0
        # fixit(x)=g(x)
        i = 2
        while i <= n
            xs[i] = try
                g(xs[i-1])
            catch
                NaN
            end
            if abs(1 - xs[i-1] / xs[i]) < ϵ || isnan(xs[i])
                break
            end
            i += 1
        end
        # gx = [r(i,x0) for i in 1:n]
        last_i = min(i, n)
        xs = filter(x -> !isnan(x), xs[1:last_i])
        ys = xs[2:end]
        xss = xs[1:end-1]
        xss, ys
    end
    function animate_fixedpoint(g, x0, ϵ; maxiters=50, limits=nothing)
        xs, ys = fixed_point(g, x0, ϵ; maxiters=50)

        plt = if isnothing(limits)
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright)
        else
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright, xlimits=limits)
        end
        anim = @animate for j ∈ 1:(length(xs)-1)
            scatter(plt, xs[1:j], repeat([0], j), label=L"x_{%$j}=%$(xs[j])")
        end
        # # annotation=[(2,5,L"x_{%$i}=%$(xs[i])",10)]
        # # gif(anim, "anim_fps15.gif", fps = 2)
        anim
    end

end

# ╔═╡ f35de386-da3a-4cbf-89ae-3049218531df
let
    # anim = animate_fixedpoint(x->x^2-2,1.1,0.001)
    # gif(anim, "anim1_fps15.gif", fps = 2)
    # g(x)=x^2-2;
    # plot([g,x->x],framestyle=:origin, label=[L"y=x^2-2"  L"y=x"])
    # scatter!([-1,2],g.([-1,2]), label="fixed points")
end

# ╔═╡ 22cb99e8-a5e4-4a16-8670-1ef1ef6f7b39
let
    g(x) = (x^2 - 1) / 3
    plot([g, x -> x], framestyle=:origin, label=[L"y=%$(g(x))" L"y=x"])
    u1, = solveit(x -> g(x) - x, -2)
    u2, = solveit(x -> g(x) - x, 2)

    scatter!([u1, u2], g.([u1, u2]), label="fixed points")
end

# ╔═╡ 9ac51023-96ca-4304-a8b0-af36c3c8f60e
md"## Fixed-Point Iteration"

# ╔═╡ 8cf9969b-c6ce-46af-9770-effd72bdf06c
cm"""
```math
\begin{array}{lcl}
g_1(x)&=&x-x^3-4x^2+10\\
	g_2(x)&=&\displaystyle \sqrt{\frac{10}{x}-4x}\\
	g_3(x)&=&\displaystyle \frac{1}{2}\sqrt{10-x^3}\\
	g_4(x)&=&\displaystyle \sqrt{\frac{10}{4+x}}\\
	g_5(x)&=&\displaystyle x-\frac{x^3+4x^2-10}{3x^2+8x}\\
\end{array}
"""

# ╔═╡ 4ed1a6ca-0d21-4257-87db-f39fff0d208f
let
    p0 = 1.5
    n = 30
    ϵ = 1e-9
    g1(x) = x - x^3 - 4x^2 + 10
    g2(x) = sqrt((10 / x) - 4x)
    g3(x) = (1 / 2.0) * sqrt(10 - x^3)
    g4(x) = sqrt(10 / (4 + x))
    g5(x) = x - (x^3 + 4x^2 - 10) / (3x^2 + 8x)
    xs1, ys1 = fixed_point(g1, p0, ϵ; maxiters=n)
    xs2, ys2 = fixed_point(g2, p0, ϵ; maxiters=n)
    xs3, ys3 = fixed_point(g3, p0, ϵ; maxiters=n)
    xs4, ys4 = fixed_point(g4, p0, ϵ; maxiters=n)
    xs5, ys5 = fixed_point(g5, p0, ϵ; maxiters=n)

    T = hcat(0:n,
        vcat(xs1, repeat([" "], n + 1 - length(xs1))),
        vcat(xs2, repeat([" "], n + 1 - length(xs2))),
        vcat(xs3, repeat([" "], n + 1 - length(xs3))),
        vcat(xs4, repeat([" "], n + 1 - length(xs4))),
        vcat(xs5, repeat([" "], n + 1 - length(xs5))))
    # T = hcat(reshape(repeat([" "],31),31,1))
    # length(xs1)
    pretty_table(HTML, T, header=vcat(:n, [Symbol("g$i") for i in 1:5]...))
    # anim2 = animate_fixedpoint(g3,p0,ϵ;maxiters=n,limits=(0,2))
    # gif(anim2, "anim2_fps15.gif", fps = 2)
end

# ╔═╡ 54b12ace-4743-49a5-9e43-8580ee43ca6a
cm"""
- __Question__: How can we find a fixed-point problem that produces a sequence that reliably and rapidly converges to a solution to a given root-finding problem?
"""



# ╔═╡ 7ecbc555-7c10-4002-a662-b3de16611269
cm"""
- __Answer__: Manipulate the root-finding problem into a fixed point problem that satisfies the conditions of Fixed-Point Theorem 2.4 and has a derivative that is as small as possible near the fixed point.
"""

# ╔═╡ 9d6e52a1-d573-4b50-a584-5e517991e8ed
md"""# 2.3 Newton's Method and its Extensions 
__(Newton’s method and Secant method only)__
"""

# ╔═╡ e8756fef-22eb-48d9-b430-24b2af6dea4f
begin
    function newton_method(f, f_prime, x0; tol=1e-7, max_iter=1000)
        x = x0
        T = Vector{Number}(undef, max_iter)

        for i in 1:max_iter
            fx = f(x)
            fpx = f_prime(x)
            T[i] = x
            # Check if the derivative is zero to avoid division by zero
            if abs(fpx) < tol
                # println("Derivative too small; stopping iteration to avoid division by zero.")
                return x, T[1:i], :small_derivative
            end

            x_new = x - fx / fpx

            # Check for convergence
            if abs(x_new - x) < tol
                T[i+1] = x_new
                return x_new, T[1:i+1], :converged
            end

            x = x_new
        end

        # println("Maximum number of iterations reached without convergence.")
        return x, T, :maxiters_reached
    end
    function animate_newton(g, gp, x0; tol=1e-7, max_iter=1000, limits=nothing)
        _, T, _ = newton_method(g, gp, x0; tol=tol, max_iter=max_iter)
        xs = T
        plt = if isnothing(limits)
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright)
        else
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright, xlimits=limits)
        end
        anim = @animate for j ∈ 1:length(xs)
            scatter(plt, xs[1:j], repeat([0], j), label=L"x_{%$j}=%$(xs[j])")
        end
        # # annotation=[(2,5,L"x_{%$i}=%$(xs[i])",10)]
        # # gif(anim, "anim_fps15.gif", fps = 2)
        anim
    end
end


# ╔═╡ d70b02b7-9697-428f-8eb7-e75f329da362
md"## Newton's Method (Newton-Raphson)"

# ╔═╡ ffb5767a-24fe-411f-abaa-baac67eaaa3d
cm"""
__Stopping Criteria__ Same as Bisection method. Nameley, select a tolerance ``\varepsilon>0`` and construct ``p_1, \ldots p_N`` until
```math
\begin{aligned}
& \left|p_N-p_{N-1}\right|<\varepsilon \\
& \frac{\left|p_N-p_{N-1}\right|}{\left|p_N\right|}<\varepsilon, \quad p_N \neq 0
\end{aligned}
```
or
```math
\left|f\left(p_N\right)\right|<\varepsilon .
```
"""

# ╔═╡ e4f0570b-0913-4f40-8a9b-8afb5ce7cbcd
let
    # x0=0
    # f(x)=cos(x)-x
    # fp(x)=-sin(x)-1
    # xsf,ysf=fixed_point(x->cos(x),x0,1e-6;maxiters=1000)
    # anim_fixed = animate_fixedpoint(x->cos(x),x0,1e-6;maxiters=1000)
    # gif(anim_fixed,"anim_fixed_n1.gif",fps=2)
    # x,xsn,flag = newton_method(f,fp,x0;tol=1e-6,max_iter=1000)
    # anim_newton = animate_newton(f,fp,x0;tol=1e-6,max_iter=1000,limits=(0.5,1))
    # gif(anim_newton,"anim_newton_n1.gif",fps=2)
    # common_len=max(length(xsf),length(xsn))
    # T = Matrix{Union{Missing,Number,String}}(missing,common_len,3)
    # T[1:common_len,1]=collect(1:common_len)
    # T[1:length(xsf),2]=xsf
    # T[1:length(xsn),3]=xsn
    # T[:,2] = map(x->ismissing(x) ?  " " : x,T[:,2])
    # T[:,3] = map(x->ismissing(x) ?  " " : x,T[:,3])
    # pretty_table(HTML,T;header=["n" ,"Fixed Point", "Newton"])
end

# ╔═╡ 8ce9ee8c-cca0-4ff5-a5a0-14991987feb0
md"## The Secant Method"

# ╔═╡ e0ab7d80-17fc-478e-90a1-4f0922bfd728
begin
    function secant_method(f, p0, p1, TOL, N0)
        i = 1
        q0 = f(p0)
        q1 = f(p1)
        T = Matrix{Number}(undef, N0 + 1, 3)  # Initialize matrix T with NaN values
        FLAG = 0  # Initialize FLAG as failure
        T[1, :] = [0, p0, q0]
        T[2, :] = [1, p1, q1]
        while i <= N0
            p = p1 - q1 * (p1 - p0) / (q1 - q0)
            fp = f(p)
            T[i+1, :] = [i, p, fp]

            if abs(p - p1) < TOL
                TT = T[1:i+1, :]  # Return the matrix up to the current iteration
                FLAG = 1  # Set FLAG as success
                return p, TT, FLAG
            end

            i += 1
            p0, q0 = p1, q1
            p1, q1 = p, fp
        end

        TT = T  # Return the matrix up to the last iteration
        p = nothing  # Indicate no valid solution found
        return p, TT, FLAG
    end
    function animate_secant(g, x0, x1, tol, max_iter; limits=nothing)
        _, T, _ = secant_method(g, x0, x1, tol, max_iter)
        xs = T[:, 2]

        plt = if isnothing(limits)
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright)
        else
            plot(g, label=L"g(x)", framestyle=:origin, legend=:topright, xlimits=limits)
        end
        anim = @animate for j ∈ 1:length(xs)
            scatter(plt, xs[1:j], repeat([0], j), label=L"x_{%$j}=%$(xs[j])")
        end
        # # annotation=[(2,5,L"x_{%$i}=%$(xs[i])",10)]
        # # gif(anim, "anim_fps15.gif", fps = 2)
        anim
    end

end

# ╔═╡ f91c04c2-4e85-4709-a76e-639929d54abd
let
    x, TT, flag = secant_method(x -> cos(x) - x, 0, 1, 1e-8, 15)
    sccess = [i == size(TT, 1) ? flag == 1 ? "Success" : "Fail" : "" for i in 1:size(TT, 1)]
    pretty_table(HTML, hcat(TT, sccess), header=["n", "p_n", "f(p_n)", "FLAG"])
end

# ╔═╡ 65193530-df0b-46f5-8653-5379a16ab779
# let
# 	anim_sec = animate_secant(x->cos(x)-x,0,1,1e-8,30;limits=(0.5,1.0))
# 	gif(anim_sec, "anim_sec_fps15.gif", fps = 2)
# end

# ╔═╡ c638088a-c6fe-406e-8ef5-f3511319aef9
md"# 3.1 Interpolation and the Lagrange Polynomials"

# ╔═╡ 70ad83fd-6df4-47b7-b5c0-0225c552d2e7
md"## Lagrange Interpolating Polynomials"

# ╔═╡ 261a5231-b450-4945-ba85-28c737e5f46f
begin
    L(y, i) = x -> prod(i == j ? 1 : (x - y[j]) / (y[i] - y[j]) for j in 1:length(y))
    LagrangeP(xs, ys) = x -> sum(ys[i] * L(xs, i)(x) for i in 1:length(xs))
    LagrangeR(xs, ndf) = (ζ, x) -> ndf(ζ) * prod((x - xs[i]) for i in 1:length(xs)) / (factorial(length(xs)))
end

# ╔═╡ 3e3022c4-06ad-4878-aa53-79e274dd40ec
cm"""

Now, for tow given point ``\left(x_0, y_0\right)`` and ``\left(x_1, y_1\right)``, define the functions
```math
L_0(x)=\frac{x-x_1}{x_0-x_1} \quad \text { and } \quad L_1(x)=\frac{x-x_0}{x_1-x_0} .
```

The __linear Lagrange interpolating polynomial__ through ``\left(x_0, y_0\right)`` and ``\left(x_1, y_1\right)`` is
```math
P(x)=L_0(x) f\left(x_0\right)+L_1(x) f\left(x_1\right)=\frac{x-x_1}{x_0-x_1} f\left(x_0\right)+\frac{x-x_0}{x_1-x_0} f\left(x_1\right) .
```
"""

# ╔═╡ 972e9347-72c0-4a53-b20a-205086a29f8b
let
    points = [(2, 4), (5, 1)]
    xs = first.(points)
    ys = last.(points)
    P(x) = LagrangeP(xs, ys)(x)
    expand(P(x))
    plt = plot(P)
    scatter(plt, xs, ys)
end

# ╔═╡ 123cec56-4f6d-445d-9f53-629161a1d487
cm"""
__General Case__:

To generalize the concept of linear interpolation, consider the construction of a polynomial of degree at most ``n`` that passes through the ``n+1`` points
```math
\left(x_0, f\left(x_0\right)\right),\left(x_1, f\left(x_1\right)\right), \ldots,\left(x_n, f\left(x_n\right)\right) \text {. }
```
(See Figure 3.4.)
"""

# ╔═╡ 26847ae3-9b6d-46b4-84e0-225177041c15
cm"""
In this case, we first construct, for each ``k=0,1, \ldots, n``, a function ``L_{n, k}(x)`` with the property that 
- ``L_{n, k}\left(x_i\right)=0`` when ``i \neq k`` and 
- ``L_{n, k}\left(x_k\right)=1``. 

To satisfy ``L_{n, k}\left(x_i\right)=0`` for each ``i \neq k`` requires that the numerator of ``L_{n, k}(x)`` contain the term
```math
\left(x-x_0\right)\left(x-x_1\right) \cdots\left(x-x_{k-1}\right)\left(x-x_{k+1}\right) \cdots\left(x-x_n\right) .
```

To satisfy ``L_{n, k}\left(x_k\right)=1``, the denominator of ``L_{n, k}(x)`` must be this same term but evaluated at ``x=x_k``. Thus,
```math
L_{n, k}(x)=\frac{\left(x-x_0\right) \cdots\left(x-x_{k-1}\right)\left(x-x_{k+1}\right) \cdots\left(x-x_n\right)}{\left(x_k-x_0\right) \cdots\left(x_k-x_{k-1}\right)\left(x_k-x_{k+1}\right) \cdots\left(x_k-x_n\right)} =\prod_{\substack{i=0 \\ i \neq k}}^n \frac{\left(x-x_i\right)}{\left(x_k-x_i\right)}.
```

A sketch of the graph of a typical ``L_{n, k}`` (when ``n`` is even) is shown in Figure 3.5.
"""

# ╔═╡ 06126200-59c5-4277-bb0c-f156ed90f51f
cm"""
The polynomial  given by
```math
P(x)=f\left(x_0\right) L_{n, 0}(x)+\cdots+f\left(x_n\right) L_{n, n}(x)=\sum_{k=0}^n f\left(x_k\right) L_{n, k}(x),
```
is called the ``\boldsymbol{n}``th __Lagrange interpolating polynomial__.
-  It is unique.
- From now on we write ``L_{k}(x)`` instead of ``L_{n, k}(x)``.
"""

# ╔═╡ c101ca1b-13af-4d0b-806b-14481f81b13e
let
    f(x) = 1 / x
    xs = [2; 2.75; 4]
    ys = f.(xs)
    P(x) = LagrangeP(xs, ys)(x) |> expand
    # expand(P(x))

    # P(4)
    plot([f, P], framestyle=:zeros, xlimit=(1, 10), label=[L"%$(f(x))" L"%$(P(x))"])
end

# ╔═╡ 4f9efa7b-a9ea-4012-87e9-7d0cedb9be54
let
    @syms ζ::Real
    f(x) = 1 / x
    xs = [2; 2.75; 4]
    ys = f.(xs)
    P(x) = LagrangeP(xs, ys)(x) |> expand
    ndf(x) = -6 * x^(-4)
    R(ζ) = LagrangeR(xs, ndf)(ζ, x)
    R(ζ)
    # P(3)
    # plot([f,P],framestyle=:zeros, xlimit=(2,4), label=[L"%$(f(x))" L"%$(P(x))"])
end

# ╔═╡ e8f8c6c1-184c-4edf-a0b7-ce13a6eaf2a0
md"# 3.3 Divided Differences"

# ╔═╡ cfbdaefd-bdcc-43e5-bdb4-87b37071c3c8
cm"""
- Suppose that ``P_n(x)`` is the ``n``th interpolating polynomial that agrees with the function ``f`` at the distinct numbers ``x_0, x_1, \ldots, x_n``. 

- The __divided differences__ of ``f`` with respect to ``x_0, x_1, \ldots, x_n`` are used to express ``P_n(x)`` in the form
```math
P_n(x)=a_0+a_1\left(x-x_0\right)+a_2\left(x-x_0\right)\left(x-x_1\right)+\cdots+a_n\left(x-x_0\right) \cdots\left(x-x_{n-1}\right),
```
``\text{-- }``for appropriate constants ``a_0, a_1, \ldots, a_n``. 
- To determine the first of these constants, ``a_0``, note that if ``P_n(x)`` is written in the form of Eq. (3.5), then evaluating ``P_n(x)`` at ``x_0`` leaves only the constant term ``a_0``; that is,
```math
a_0=P_n\left(x_0\right)=f\left(x_0\right) \text {. }
```
- Similarly, when ``P(x)`` is evaluated at ``x_1``, the only nonzero terms in the evaluation of ``P_n\left(x_1\right)`` are the constant and linear terms,
```math
f\left(x_0\right)+a_1\left(x_1-x_0\right)=P_n\left(x_1\right)=f\left(x_1\right) ;
```

"""

# ╔═╡ e3deaaf2-9c9b-4693-b903-4c3ecfda1cd7
cm"""
so,
```math
a_1=\frac{f\left(x_1\right)-f\left(x_0\right)}{x_1-x_0} .
```

> We now introduce the __divided-difference notation__, which is related to Aitken's __``\Delta^2`` notation__. 

- The __zeroth divided difference__ of the function ``f`` with respect
to ``x_i``, denoted ``f\left[x_i\right]``, is simply the value of ``f`` at ``x_i`` :
```math
f\left[x_i\right]=f\left(x_i\right) .
```

- the __first divided difference__ of ``f`` with respect to ``x_i`` and ``x_{i+1}`` is denoted ``f\left[x_i, x_{i+1}\right]`` and defined as
```math
f\left[x_i, x_{i+1}\right]=\frac{f\left[x_{i+1}\right]-f\left[x_i\right]}{x_{i+1}-x_i} .
```

- The __second divided difference__, ``f\left[x_i, x_{i+1}, x_{i+2}\right]``, is defined as
```math
f\left[x_i, x_{i+1}, x_{i+2}\right]=\frac{f\left[x_{i+1}, x_{i+2}\right]-f\left[x_i, x_{i+1}\right]}{x_{i+2}-x_i} .
```
-  Similarly, after the ``(k-1)`` st divided differences,
```math
f\left[x_i, x_{i+1}, x_{i+2}, \ldots, x_{i+k-1}\right] \text { and } f\left[x_{i+1}, x_{i+2}, \ldots, x_{i+k-1}, x_{i+k}\right] \text {, }
```
have been determined, the ``k`` th divided difference relative to ``x_i, x_{i+1}, x_{i+2}, \ldots, x_{i+k}`` is
```math
f\left[x_i, x_{i+1}, \ldots, x_{i+k-1}, x_{i+k}\right]=\frac{f\left[x_{i+1}, x_{i+2}, \ldots, x_{i+k}\right]-f\left[x_i, x_{i+1}, \ldots, x_{i+k-1}\right]}{x_{i+k}-x_i}
```
"""

# ╔═╡ e68ab39e-fdb5-4f39-8f18-c25f65e99c05
cm"""
- The process ends with the single nth divided difference,
```math
f\left[x_0, x_1, \ldots, x_n\right]=\frac{f\left[x_1, x_2, \ldots, x_n\right]-f\left[x_0, x_1, \ldots, x_{n-1}\right]}{x_n-x_0}
```
"""

# ╔═╡ 3d5e1274-386b-4511-8353-0188b2b65eb0
begin
    function newton_devided_diff(x, y)
        n = length(x)
        F = Matrix{Real}(undef, n, n)
        F[:, 1] = y
        for col in 2:n
            for row in col:n
                F[row, col] = (F[row, col-1] - F[row-1, col-1]) / (x[row] - x[row-col+1])
            end
        end
        diag(F)
        # F
    end
    function NDDP(xs, ys)
        F = newton_devided_diff(xs, ys)
        n = length(xs)
        return x -> begin
            F[1] + sum(F[i] * reduce(*, map(u -> (x - u), xs[1:i-1])) for i in 2:n)
        end
    end
end

# ╔═╡ 5d5adcae-2ad4-44b7-af71-a8ee276df8a3
let
    xs = [1.0; 1.3; 1.6; 1.9; 2.2]
    ys = [0.7651977; 0.6200860; 0.4554022; 0.2818186; 0.1103623]
    F = newton_devided_diff(xs, ys)
    P1(x) = NDDP(xs, ys)(x) |> expand
    P1(x)
    P2(x) = LagrangeP(xs, ys)(x) |> expand
    Dict(
        :newton_devided_diff => P1(x),
        :lagrange => P2(x)
    )
    # expand(P(x))
    P1(1.5), P2(1.5)
end

# ╔═╡ c3fb95a4-c9c0-4fc6-85d4-f6cf8a8258e4
md"""
# 3.5 Cubic Spline Interpolation 
## Piecewise-Polynomial Approximation
"""

# ╔═╡ 61b0d517-94fc-4c6c-814b-8dfe51da6f1b
md"## Cubic Splines"

# ╔═╡ 5f768b64-c67e-4470-bf44-35bb539e1441
begin
    function natural_spline(x, y)
        n = length(x)
        h = map(i -> x[i+1] - x[i], 1:n-1)
        a = copy(y)
        α = map(2:n-1) do i
            (3 / h[i]) * (a[i+1] - a[i]) - (3 / h[i-1]) * (a[i] - a[i-1])
        end
        l = Vector{Real}(undef, n)
        μ = Vector{Real}(undef, n - 1)
        z = copy(l)
        c = Vector{Real}(undef, n)
        b = copy(μ)
        d = copy(μ)
        l[1] = 0
        μ[1] = 0
        z[1] = 0
        l[n] = 1
        z[n] = 0
        foreach(2:n-1) do i
            l[i] = 2(x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
            μ[i] = h[i] / l[i]
            z[i] = (α[i-1] - h[i-1] * z[i-1]) / l[i]
        end

        c[n] = 0
        foreach(n-1:-1:1) do i
            c[i] = z[i] - μ[i] * c[i+1]
            b[i] = (a[i+1] - a[i]) / h[i] - (h[i] * (c[i+1] + 2c[i]) / 3)
            d[i] = (c[i+1] - c[i]) / (3h[i])
        end
        a, b, c, d
    end
end

# ╔═╡ c4a80b01-c1e3-4196-8408-3054d0aca71e
md"## Construction of a Cubic Spline"

# ╔═╡ 25317d39-d107-41e4-a5e3-5ee11a0f8bd5
cm"""
As the preceding example demonstrates, 
- a spline defined on an interval that is divided into ``n`` subintervals will require determining __``4 n``__ constants. To construct the cubic spline interpolant for a given function ``f``, the conditions in the definition are applied to the cubic polynomials
```math
S_j(x)=a_j+b_j\left(x-x_j\right)+c_j\left(x-x_j\right)^2+d_j\left(x-x_j\right)^3,
\quad j \in \{0,1,2,\cdots,n-1\} 
```
"""

# ╔═╡ d65a3a75-9e63-4087-912a-89e5e973a171
md"## Natural Splines"

# ╔═╡ bff8bd64-7a8c-4105-96b1-a3f8e8f0cd41
cm"""
- Natural Cubic Spline produces a linear system
```math
A x = b
```
where 
```math
A = \begin{bmatrix}
1 & 0 & 0 & \cdots &\cdots & 0 & 0 \\
h_0 & 2(h_0 + h_1) & h_1 & \cdots &\cdots & 0 & 0 \\
0 & h_1 & 2(h_1 + h_2) & \cdots &\cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \cdots &\vdots & \vdots \\
0 & 0 & 0 & h_{n-1} & \cdots & 2(h_{n-2} + h_{n-1}) & h_{n-1} \\
0 & 0 & 0 & \cdots &\cdots &  & 1 \\
\end{bmatrix}
```

and the vectors ``x`` and ``b`` are defined as

```math
\mathbf{b}=\left[\begin{array}{c}0 \\ \frac{3}{h_1}\left(a_2-a_1\right)-\frac{3}{h_0}\left(a_1-a_0\right) \\ \vdots \\ \frac{3}{h_{n-1}}\left(a_n-a_{n-1}\right)-\frac{3}{h_{n-2}}\left(a_{n-1}-a_{n-2}\right) \\ 0\end{array}\right] \quad \text{ and} \quad \mathbf{x}=\left[\begin{array}{c}c_0 \\ c_1 \\ \vdots \\ c_n\end{array}\right]
```
"""

# ╔═╡ c5f4f53e-6200-4a02-ada7-72e8ae60e38b
cm"""
- ``a_j = f(x_j)`` for all ``j=0,\cdots, n``
- ``b_j=\displaystyle\frac{1}{h_j}\left(a_{j+1}-a_j\right)-\frac{h_j}{3}\left(2 c_j+c_{j+1}\right)`` for all ``j=0,\cdots, n-1``
- ``d_j = \displaystyle \frac{c_{j+1}-c_j}{3h_j}`` for all ``j=0,\cdots, n-1``

"""

# ╔═╡ c2c69250-1f18-4c6d-b8b7-3de26ba0818a
let
    xs = [1, 2, 3]
    ys = [2, 3, 5]
    a, b, c, d = natural_spline(xs, ys)
    Si(i) = x -> a[i] + b[i] * (x - xs[i]) + c[i] * (x - xs[i])^2 + d[i] * (x - xs[i])^3
    S(x) = begin
        i = findfirst(j -> xs[j] <= x <= xs[j+1], 1:length(xs)-1)
        Si(i)(x)
    end
    Si(2)(x)

    # A = [1 0 0;1 4 0;0 0 1]
    # b = [0;3;0]
    # c = A\b
    # h = [1;1;1]
    # a = ys
    # n = length(xs)
    # b = [(a[i+1]-a[i])/(h[i]) - h[i]*(c[i+1]+2c[i])/3 for i in 1:n-1]
    # d = [(c[i+1]-c[i])/(3h[i]) for i in 1:n-1]
    # X = [x;x]
    # S = a[1:2] .+ b.*(X.-xs[1:2]) .+  c[1:2] .* (X.-xs[1:2]).^2 .+  d .* (X.-xs[1:2]).^3

end

# ╔═╡ e45f6f80-7a93-4064-97b0-d6e7269a1b68
let
    xs = 0:3
    ys = exp.(xs)
    a, b, c, d = natural_spline(xs, ys)

    Si(i) = x -> a[i] + b[i] * (x - xs[i]) + c[i] * (x - xs[i])^2 + d[i] * (x - xs[i])^3
    S(x) = begin
        i = findfirst(j -> xs[j] <= x <= xs[j+1], 1:length(xs)-1)
        Si(i)(x)
    end
    plot(0:0.0001:3, [exp.(0:0.0001:3), S.(0:0.0001:3)], label=[L"f(x)=e^x" L"S(x)"])

    a = copy(ys)
    n = length(xs)
    h = [xs[i+1]-xs[i] for i in 1:n-1]
    Aii = map(i-> i in (1,n) ? 1 : 2*(h[i-1]+h[i]), 1:n)
    Ait = map(i-> i==1 ? 0 : h[i], 1:n-1)
    Aib = map(i-> i==(n-1) ? 0 : h[i], 1:n-1)
    A = diagm(0=>Aii,1=>Ait,-1=>Aib)
    ho3 = 3 ./ h
    rhs = map(i-> i in (1,n) ? 0 : ho3[i]*(a[i+1]-a[i])-ho3[i-1]*(a[i]-a[i-1]),1:n)
    c = A\rhs

    b = [(a[i+1]-a[i])/(h[i]) - h[i]*(c[i+1]+2c[i])/3 for i in 1:n-1]
    d = [(c[i+1]-c[i])/(3h[i]) for i in 1:n-1]
    X = repeat([x],n-1)
    S = a[1:n-1] .+ b.*(X.-xs[1:n-1]) .+  c[1:n-1] .* (X.-xs[1:n-1]).^2 .+  d .* (X.-xs[1:n-1]).^3
end

# ╔═╡ 3115c2f6-8574-4d09-af83-db8e09731d07
L"""
\int_0^3 e^ x dx = \left. e^x \right|_0^3 = e^3 - 1 = %$(exp(3)-1)
"""

# ╔═╡ 01b33b12-1dad-45e9-a191-8f9619b2fac5
L"""
\int_0^3 e^ x dx \approx \int_0^1 S_0(x) dx  +\int_1^2 S_1(x) dx+\int_2^3 S_2(x) dx
"""

# ╔═╡ b39bcd45-dfc3-4ddb-95e8-7272ba591682
let
    xs = 0:3
    ys = exp.(xs)
    a, b, c, d = natural_spline(xs, ys)
    Si(i) = x -> a[i] + b[i] * (x - xs[i]) + c[i] * (x - xs[i])^2 + d[i] * (x - xs[i])^3
    integratit(i, a, b) = solve(IntegralProblem((x, _) -> Si(i)(x), (a, b)), QuadGKJL()).u
    L"""
    \int_0^3 e^ x dx \approx %$(integratit(1,0,1))  +%$(+integratit(2,1,2))+%$(+integratit(3,2,3)) = %$(round((integratit(1,0,1)+integratit(2,1,2)+integratit(3,2,3)),digits=4))
    """

end

# ╔═╡ eb60c5be-bad0-4c39-9842-3210dbaba8d3
md"## Clamped Splines"

# ╔═╡ 0c695eed-ebe3-4a26-9a32-7054bf1a51c3
let
    xs = [1; 2; 3]
    ys = [2; 3; 5]
    A = [
        1 1 1 0 0 0
        0 0 0 1 1 1
        1 2 3 -1 0 0
        0 1 3 0 -1 0
        0 1 0 0 0 0
        0 0 0 0 1 3
    ]
    b = [
        1
        2
        0
        0
        0
        0
    ]
    b0, c0, d0, b1, c1, d1 = A \ b
    A2 = [
        1 1 1 0 0 0
        0 0 0 1 1 1
        1 2 3 -1 0 0
        0 1 3 0 -1 0
        1 0 0 0 0 0
        0 0 0 1 2 3
    ]
    b2 = [
        1
        2
        0
        0
        2
        1
    ]
    b01, c01, d01, b11, c11, d11 = A2 \ b2
    # A2\b2
    aN = [2; 3; 5]
    bN = [b0; b1]
    cN = [c0; c1]
    dN = [d0; d1]
    SiN(i) = x -> aN[i] + bN[i] * (x - xs[i]) + cN[i] * (x - xs[i])^2 + dN[i] * (x - xs[i])^3
    S(x) = begin
        i = findfirst(j -> xs[j] <= x <= xs[j+1], 1:length(xs)-1)
        SiN(i)(x)
    end
    bC = [b01; b11]
    cC = [c01; c11]
    dC = [d01; d11]
    SiC(i) = x -> aN[i] + bC[i] * (x - xs[i]) + cC[i] * (x - xs[i])^2 + dC[i] * (x - xs[i])^3
    SC(x) = begin
        i = findfirst(j -> xs[j] <= x <= xs[j+1], 1:length(xs)-1)
        SiC(i)(x)
    end
    plot(1:0.01:3, [x -> S(x), x -> SC(x)])
	map(f->f(x),SiC.(1:2))
end

# ╔═╡ 92244fae-7546-4608-83ad-6dfe7425243c
let
    A = [
        1 1 1 0 0 0
        0 0 0 1 1 1
        1 2 3 -1 0 0
        0 1 3 0 -1 0
        1 0 0 0 0 0
        0 0 0 1 2 3
    ]
    b = [
        1
        2
        0
        0
        2
        1
    ]
    b0, c0, d0, b1, c1, d1 = A \ b

end

# ╔═╡ 5fff03aa-5063-4bd0-bad1-6a3e568b25a1
md"# 8.1 Discrete Least Squares Approximation"

# ╔═╡ f50ef025-5faa-4a17-9307-4687bd3daac6
begin
	struct FitModel
		name::Symbol
		n::Union{Int,Nothing}
	end
	function discrete_least_squares(xs,ys;model::FitModel=FitModel(:poly,1))
		m = length(xs)
		uxs = unique(xs)
		if length(uxs)!=length(xs)
			@error "The x data must be unique"
			return nothing
		end
		fit_model = model.name
		data = if fit_model==:power
			if any(ys.<0) || any(xs.<0)	
				return nothing, "The data must be potivive", nothing
			end
			log.(xs),log.(ys),1
		elseif fit_model== :exp
			if any(ys.<0) 
				return nothing, "The data must be potivive", nothing
			end
			xs,log.(ys),1
		else
			nn = isnothing(model.n) ? 1 : model.n
			xs,ys, nn 
		end
		nxs,nys, n = data
		if isnothing(nxs)
			@error nys
			return
		end
		sz = n + 1
		entries = map(0:2n) do i
			sum(nxs.^i)
		end
		
		A =stack(map(1:sz) do i
			entries[i:i-1+sz]
		end)
		
		
		rhs = map(0:n) do i
			sum(nys.*nxs.^i)
		end
		A\rhs
	end	
end

# ╔═╡ 4fde53d9-ae78-495e-9898-65af518d6b7f
cm"""
Consider the problem of estimating the values of a function at nontabulated points, given the experimental data in Table 
| ``x_i`` | ``y_i`` |
| ---: | ---: |
| 1 | 1.3 |
| 2 | 3.5 |
| 3 | 4.2 |
| 4 | 5.0 |
| 5 | 7.0 |
| 6 | 8.8 |
| 7 | 10.1 |
| 8 | 12.5 |
| 9 | 13.0 |
| 10 | 15.6 |

"""

# ╔═╡ 44e027ce-9f36-4ffb-bf81-b10905107770
cm" 1. __Let's plot these points__"

# ╔═╡ 0cfec2e7-9166-45ad-b7cd-9585a1205e22
let
    xs = 1:10
    ys = [1.3, 3.5, 4.2, 5.0, 7.0, 8.8, 10.1, 12.5, 13.0, 15.6]
    scatter(xs, ys, label=nothing)
end

# ╔═╡ 2c65f0af-38dd-4d5b-868d-c62b6f2653f9
cm"""
2. __What is the relationship between ``x`` and ``y``?__ (Linear)
3. We need to find ``f(x) = a_1 x + a_0`` that approximates these data and gives a good predictor.
"""

# ╔═╡ 9219bf2c-6e18-4a82-bf59-5050701bc404
begin
    a0Slider = @bind a0_slider NumberField(-2:0.001:15, default=0)
    a1Slider = @bind a1_slider NumberField(0:0.001:5, default=1)
    cm"""
    ``a_0`` = $(a0Slider)

    ``a_1`` = $(a1Slider)

    """
end

# ╔═╡ f24491c0-11a1-45b0-adac-cb9ce23f300f
let
    least_squared(ys, fxs) = round(sum((ys - fxs) .^ 2), digits=4)
    xs = 1:10
    ys = [1.3, 3.5, 4.2, 5.0, 7.0, 8.8, 10.1, 12.5, 13.0, 15.6]
    p1 = scatter(xs, ys, label=nothing, framestyle=:origin, xlimits=(0, 15), ylimits=(-1, 20))
    least_line(x) = a0_slider + a1_slider * x
    p1 = plot(p1, x -> least_line(x), label=L"y=%$(a0_slider)+%$(a1_slider)*x")
    scatter(p1, xs, x -> least_line(x), label=L"points on line")
    annotate!([(2.2, 15, "Error= " * L"%$(least_squared(least_line.(xs),ys))")])
end

# ╔═╡ 941815a4-e9e9-4d9d-b019-31da9ff3b176
let
	xs = 1:10
	ys = [1.3, 3.5, 4.2, 5.0, 7.0, 8.8, 10.1, 12.5, 13.0, 15.6]
	xs_sum = sum(xs)
	ys_sum = sum(ys)
	xs2_sum = sum(xs .^ 2)
	xsys_sum = sum(xs .* ys)
	A = [length(xs)+1 xs_sum
		xs_sum xs2_sum
	]
	b = [ys_sum; xsys_sum]
	a₀,a₁ = A \ b
end

# ╔═╡ 35608f80-9cf7-441b-a4f8-b5617899528d
md"## Polynomial Least Squares"

# ╔═╡ 4b236b75-12dd-4638-b968-4923123af93d
cm"""
The general problem of approximating a set of data, 
``\left\{\left(x_i, y_i\right) \mid i=1,2, \ldots, m\right\}``, with an algebraic polynomial

```math
P_n(x)=a_n x^n+a_{n-1} x^{n-1}+\cdots+a_1 x+a_0,
```

of degree ``n < m-1``, using the least squares procedure is handled similarly. We choose the constants ``a_0, a_1, \ldots, a_n`` to minimize the least squares error ``E=E_2\left(a_0, a_1, \ldots, a_n\right)``, where
```math
\begin{aligned}
E & =\sum_{i=1}^m\left(y_i-P_n\left(x_i\right)\right)^2 \\
& =\sum_{i=1}^m y_i^2-2 \sum_{i=1}^m P_n\left(x_i\right) y_i+\sum_{i=1}^m\left(P_n\left(x_i\right)\right)^2 \\
& =\sum_{i=1}^m y_i^2-2 \sum_{i=1}^m\left(\sum_{j=0}^n a_j x_i^j\right) y_i+\sum_{i=1}^m\left(\sum_{j=0}^n a_j x_i^j\right)^2 \\ & =\sum_{i=1}^m y_i^2-2 \sum_{j=0}^n a_j\left(\sum_{i=1}^m y_i x_i^j\right)+\sum_{j=0}^n \sum_{k=0}^n a_j a_k\left(\sum_{i=1}^m x_i^{j+k}\right)
\end{aligned}
```

This gives ``n+1`` normal equations in the ``n+1`` unknowns ``a_j``. These are
```math
\sum_{k=0}^n a_k \sum_{i=1}^m x_i^{j+k}=\sum_{i=1}^m y_i x_i^j, \quad \text { for each } j=0,1, \ldots, n \text {. }
```

It is helpful to write the equations as follows:
```math
\begin{array}{r}
a_0 \sum_{i=1}^m x_i^0+a_1 \sum_{i=1}^m x_i^1+a_2 \sum_{i=1}^m x_i^2+\cdots+a_n \sum_{i=1}^m x_i^n=\sum_{i=1}^m y_i x_i^0, \\
a_0 \sum_{i=1}^m x_i^1+a_1 \sum_{i=1}^m x_i^2+a_2 \sum_{i=1}^m x_i^3+\cdots+a_n \sum_{i=1}^m x_i^{n+1}=\sum_{i=1}^m y_i x_i^1, \\
\vdots \\
a_0 \sum_{i=1}^m x_i^n+a_1 \sum_{i=1}^m x_i^{n+1}+a_2 \sum_{i=1}^m x_i^{n+2}+\cdots+a_n \sum_{i=1}^m x_i^{2 n}=\sum_{i=1}^m y_i x_i^n .
\end{array}
```

These normal equations have a unique solution provided that the ``x_i`` are distinct 
"""

# ╔═╡ 3141d33a-677b-47bc-ba45-fb906c984259
let
	xs = 0:0.25:1
	ys = [1.0000;1.2840;1.6487;2.1170;2.7183]
	xs_s = sum(xs)
	ys_s = sum(ys)
	xs2_s = sum(xs .^ 2)
	xs3_s = sum(xs .^ 3)
	xs4_s = sum(xs .^ 4)
	yx_s = sum(xs .* ys)
	yx2_s = sum((xs .^2) .* ys)
	A = [
		5 xs_s xs2_s 
		xs_s xs2_s xs3_s 
		xs2_s xs3_s xs4_s 
		
	]
	b =[
		ys_s
		yx_s
		yx2_s
	]
	
	a₀, a₁, a₂ = A\b
	# a₀, a₁, a₂ = discrete_least_squares(xs,ys,model=FitModel(:poly,2))
	p1 = scatter(xs,ys,mark=(:circle,:red,4), label="Data")
	plot(p1, x->a₀ + a₁*x + a₂*x^2,label=L"y=%$(a₀) + %$(a₁)x + %$(a₂)x^2",c=:purple)
end

# ╔═╡ aa7ba0e8-c1a3-4bdf-8bf5-aa389bc3795e
let
	xs = 1:0.25:2 
	ys = [5.10;5.79;6.53;7.45;8.46]
	# lnys = sum(log.(ys))
	at,b = discrete_least_squares(xs,ys,model=FitModel(:exp,nothing))
	a = exp(at)
	p1 = scatter(xs,ys,label="Data")
	plot(p1,x->a*exp(b*x),label=L"y=%$(a) e^{%$b x}")
end

# ╔═╡ 92f9bc62-9d0f-40ee-832c-dff52a5d935b
let
	xs = 1:0.25:2 
	ys = rand(-0.1:0.001:0.1,length(xs)) + 3*xs.^2.5 #[5.10;5.79;6.53;7.45;8.46]
	lnb,a = discrete_least_squares(xs,ys,model=FitModel(:power,nothing))
	b = exp(lnb)
	p1 = scatter(xs,ys,label="Data",mark=(:hexagon,7))
	plot(p1,x->b*(x^a),label=L"y=%$(b) x^{%$a}",lw=1.5)
end

# ╔═╡ 4171e8d8-5821-4a11-ba03-b8f88d785e68
md"# 4.1 Numerical Differentiation"

# ╔═╡ 479e56b4-cee6-461e-8e87-21127f8272f5
begin
	diff_2point(x0::Real,f::Function,h::Float64)=(f(x0+h)-f(x0))/h
	diff_2point(y0::Real,y1::Real,h::Float64)=(y1-y0)/h
	diff_3point(x0::Real,xj::Vector{<:Real},f::Function,h::Vector{Float64};mid::Bool=false) = begin
		map(hh->diff_3point(x0,xj,f(xj),hh;mid=mid),h)
	end
	diff_3point(x0::Real,xj::Vector{<:Real},yj::Vector{<:Real},h::Vector{Float64};mid::Bool=false) = begin
		map(hh->diff_3point(x0,xj,yj,hh;mid=mid),h)
	end
	diff_3point(x0::Real,xj::Vector{<:Real},f::Function,h::Float64;mid::Bool=false) = begin
		diff_3point(x0,xj,f(xj),h;mid=mid)
	end
	diff_3point(x0::Real,xj::Vector{<:Real},yj::Vector{<:Real},h::Float64;mid::Bool=false) = begin
		j = findfirst(x->x==x0,xj)
		sz = length(xj)
		sgn = sign(h)
		
		if mid
			k1 =findfirst(x->x==(x0+h),xj)
			k2 =findfirst(x->x==(x0-h),xj)
			if any(isnothing.([k1,k2])) 
				nothing
			else
				(1/2h)*(yj[k1]-yj[k2])				
			end
		else
			sgn=Int(sgn)
			k1 =findfirst(x->x==(x0+h),xj)
			k2 =findfirst(x->x==(x0+2h),xj)
			if any(isnothing.([k1,k2]))
				nothing
			else
				(1/2h)*(-3*yj[j]+4*yj[k1]-yj[k2])	
			end
		end
			
	end
end

# ╔═╡ e30e4b18-9d86-46d1-89f3-c15c156aff6a
cm"""
The derivative of the function ``f`` at ``x_0`` is

```math
f^{\prime}\left(x_0\right)=\lim _{h \rightarrow 0} \frac{f\left(x_0+h\right)-f\left(x_0\right)}{h} .
```

This formula gives an obvious way to generate an approximation to ``f^{\prime}\left(x_0\right)``; simply compute

```math
\frac{f\left(x_0+h\right)-f\left(x_0\right)}{h}
```
"""

# ╔═╡ 61de967c-5ba7-44db-8aad-175026d4a62d
cm"""
- To approximate ``f^{\prime}\left(x_0\right)``, suppose first that ``x_0 \in(a, b)``, where ``f \in C^2[a, b]``, and that ``x_1=x_0+h`` for some ``h \neq 0`` that is sufficiently small to ensure that ``x_1 \in[a, b]``. 
- We construct the first Lagrange polynomial ``P_{1}(x)`` for ``f`` determined by ``x_0`` and ``x_1``, with its error term:
```math
f(x)=P_{1}(x)+\frac{\left(x-x_0\right)\left(x-x_1\right)}{2!} f^{\prime \prime}(\xi(x))
```
"""

# ╔═╡ 15ba7f52-c201-4c7f-be84-da0bdf6a75f9
cm"""
- When ``x`` is ``x_0``, however, the coefficient of ``D_x f^{\prime \prime}(\xi(x))`` is 0 , and the formula simplifies to
```math
f^{\prime}\left(x_0\right)=\frac{f\left(x_0+h\right)-f\left(x_0\right)}{h}-\frac{h}{2} f^{\prime \prime}(\xi) .
```
- For small values of ``h``, the difference quotient ``\left[f\left(x_0+h\right)-f\left(x_0\right)\right] / h`` can be used to approximate ``f^{\prime}\left(x_0\right)`` with an error bounded by ``M|h| / 2``, where ``M`` is a bound on ``\left|f^{\prime \prime}(x)\right|`` for ``x`` between ``x_0`` and ``x_0+h``. 
- This formula is known as  
  - the __forward-difference formula__ if ``h>0`` and 
  - the __backward-difference formula__ if ``h<0``.

"""

# ╔═╡ 4909326c-c521-40f6-9941-ade4f56adfa7
let
	x0=1.8
	f(x)=log(x)
	h = [0.1; 0.05; 0.01;0.00000001]
	diff_2point.(x0,f,h)
	abs.(1/x0 .- diff_2point.(x0,f,h))
	# M = 1/x0^2
	# M*h/2,1/x0 .- diff_2point.(x0,f,h)
end

# ╔═╡ 9c2ec805-5efb-4b27-aa80-825220df4d66


# ╔═╡ c7f882dd-faa0-45f5-bf2e-9100c72094ba
md"## $(n+1)$-point formula to approximate $f^{\prime}\left(x_j\right)$."

# ╔═╡ 619ac568-e272-4e55-8dfc-c965f4b878c1
cm"""
To obtain general derivative approximation formulas, suppose that ``\left\{x_0, x_1, \ldots, x_n\right\}`` are ``(n+1)`` distinct numbers in some interval ``I`` and that ``f \in C^{n+1}(I)``. From Theorem 3.3 on page 109,
```math
f(x)=\sum_{k=0}^n f\left(x_k\right) L_k(x)+\frac{\left(x-x_0\right) \cdots\left(x-x_n\right)}{(n+1)!} f^{(n+1)}(\xi(x)),
```
"""

# ╔═╡ 4a7f2e34-fe1e-4a34-a942-a4cf29938167
md"## Three-Point Endpoint Formula"

# ╔═╡ ac71abdd-790e-403d-97f3-c88fe62b40b7
md"## Three-Point Midpoint Formula"

# ╔═╡ ce129375-fc38-46ae-a9b7-90605211620e
let
	xs = collect(1.8:0.1:2.2)
	df(x) = exp(x)*(1+x)
	fxs = [10.889365;12.703199;14.778112;17.148957;19.855030]
	
	# diff_3point(2.0,xs,x->x .* exp.(x),[0.2,-0.1],mid=true)
	# diff_3point(2.0,xs,fxs,[0.1,0.2],mid=true)
	  zip([:endpoint_postive_0_1,:endpoint_negative_0_1,:midpoint_0_1,:midpoint_0_2],abs.(vcat(diff_3point(2.0,xs,fxs,[0.1,-0.1],mid=false),diff_3point(2.0,xs,fxs,[0.1,0.2],mid=true)) .- df(2.0)) |> y-> map(d->@sprintf("%.2e",d),y) )   |> Dict 
	# # diff_3point(2.0,xs,fxs,[0.1],mid=false)
	
end

# ╔═╡ c85a6afe-cb7f-4188-8950-d216956b7e7a


# ╔═╡ e4e2c6aa-a849-4841-974e-a30e55843718
md"## Second Derivative Midpoint Formula"

# ╔═╡ 997136cf-24c5-45b4-bec1-87d28504d7b0
let
	xs = 1.8:0.01:2.2
	dff(x) = exp(x)*(2+x)
	fxs = xs.*exp.(xs)#[10.889365;12.703199;14.778112;17.148957;19.855030]
	adff(xj,h)= begin
		j = findfirst(x->x==xj,xs)
		k1 = findfirst(x->x==xj+h,xs)
		k2 = findfirst(x->x==xj-h,xs)
		
		(1/h^2)*(fxs[k1]+fxs[k2]-2*fxs[j])
		
	end
	adff(2.0,0.01)
end

# ╔═╡ 233e1c2d-2d45-4aa7-be03-6bc4dc76d716
md"# 4.3 Elements of Numerical Integration"

# ╔═╡ 02e6a419-66f8-4aec-b456-6667b629308a
begin
	# Trapezoidal rule and Simpson'
	
	function trapezoidal(xs::AbstractRange,f::Function)
		trapezoidal(collect(xs),f)
	end
	function trapezoidal(xs::Vector{<:Real},f::Function)
		trapezoidal(xs,f.(xs))
	end
	function trapezoidal(a::Real,b::Real,f::Function,n=1)
		h = (b-a)/n
		xs = a:h:b 
		trapezoidal(xs,f.(xs))
	end
	function trapezoidal(xs::AbstractRange,ys::Vector{<:Real})
		nxs = collect(xs)
		trapezoidal(nxs,ys)
	end
	function trapezoidal(xs::Vector{<:Real},ys::Vector{<:Real})
		n = length(xs)
		@assert  n ==length(ys) "length of `xs` must be equal to length of `ys`"
		@assert n>1 "The length of the grid must exceed 1 for Tranpezoidal rule. "
		XEND = ys[end]+ys[begin]
		XMID = sum(ys[i] for i in 2:n-1;init=0)
		h = (xs[end]-xs[begin])/(n-1)
		(h/2)*(XEND+2XMID)
	end
	function simpson(xs::AbstractRange,f::Function)
		
		simpson(collect(xs),f)
	end
	function simpson(xs::AbstractRange,ys::Vector{<:Real})
		nxs = collect(xs)
		simpson(nxs,ys)
	end
	function simpson(a::Real,b::Real,f::Function,n::Int=2)
		h = (b-a)/n
		nxs = a:h:b
		simpson(nxs,f.(nxs))
		
	end
	function simpson(xs::Vector{<:Real},f::Function)
		n = length(xs)
		
		nxs = Vector{Real}(undef,2n-1)
		nxs[1:2:end] = xs
		nxs[2:2:end-1] = 0.5*(xs[1:n-1] + xs[2:n])
		simpson(nxs,f.(nxs))
		
	end
	function simpson(xs::Vector{<:Real},ys::Vector{<:Real})
		n = length(xs)
		
		@assert n ==length(ys) "length of `xs` must be equal to length of `ys`"
			
		@assert n > 2 "Simpson needs at least 3 points"
			
		@assert isodd(n) "length of `y=f(x)` must be an odd number"
			
		XEVEN = sum(ys[i] for i in 3:2:n-1;init=0)
		XODD  = sum(ys[i] for i in 2:2:n-1;init=0)
		XEND  = ys[begin]+ys[end]
		h = (xs[end]-xs[begin])/(n-1)
		
		(h/3)*(XEND + 2XEVEN + 4XODD)
	end
		
end

# ╔═╡ 0663cf1a-a42e-440f-931f-8321cc7d65f7
cm"""
- where ``\xi(x)`` is in ``[a, b]`` for each ``x`` and
```math
a_i=\int_a^b L_i(x) d x, \quad \text { for each } i=0,1, \ldots, n .
```

- The quadrature formula is, therefore,
```math
\int_a^b f(x) d x \approx \sum_{i=0}^n a_i f\left(x_i\right)
```
with error given by
```math
E(f)=\frac{1}{(n+1)!} \int_a^b \prod_{i=0}^n\left(x-x_i\right) f^{(n+1)}(\xi(x)) d x
```
"""

# ╔═╡ e0c0c21a-8f1b-4b8e-8b4f-9fe8821c556b
md"## The Trapezoidal Rule"

# ╔═╡ 3119f567-fa6f-48a7-ac66-0cf630e3418f
cm"""

To derive the Trapezoidal rule for approximating ``\int_a^b f(x) d x``, let ``x_0=a, x_1=b, h=b-a`` and use the linear Lagrange polynomial:
```math
P_1(x)=\frac{\left(x-x_1\right)}{\left(x_0-x_1\right)} f\left(x_0\right)+\frac{\left(x-x_0\right)}{\left(x_1-x_0\right)} f\left(x_1\right) \text {. }
```

Then
```math
\begin{aligned}
\int_a^b f(x) d x= & \int_{x_0}^{x_1}\left[\frac{\left(x-x_1\right)}{\left(x_0-x_1\right)} f\left(x_0\right)+\frac{\left(x-x_0\right)}{\left(x_1-x_0\right)} f\left(x_1\right)\right] d x \\
& +\frac{1}{2} \int_{x_0}^{x_1} f^{\prime \prime}(\xi(x))\left(x-x_0\right)\left(x-x_1\right) d x .
\end{aligned}
```

The product ``\left(x-x_0\right)\left(x-x_1\right)`` does not change sign on ``\left[x_0, x_1\right]``, so the Weighted Mean Value Theorem for Integrals 1.13 can be applied to the error term to give, for some ``\xi`` in ``\left(x_0, x_1\right)``,
```math
\begin{aligned}
\int_{x_0}^{x_1} f^{\prime \prime}(\xi(x))\left(x-x_0\right)\left(x-x_1\right) d x & =f^{\prime \prime}(\xi) \int_{x_0}^{x_1}\left(x-x_0\right)\left(x-x_1\right) d x \\
& =f^{\prime \prime}(\xi)\left[\frac{x^3}{3}-\frac{\left(x_1+x_0\right)}{2} x^2+x_0 x_1 x\right]_{x_0}^{x_1} \\
& =-\frac{h^3}{6} f^{\prime \prime}(\xi) .
\end{aligned}
```
"""

# ╔═╡ e880c7b4-ab83-4bfc-ac7d-3d1616e3f925
md"## Simpson's Rule
Simpson's rule results from integrating over $[a, b]$ the second Lagrange polynomial with equally spaced nodes $x_0=a, x_2=b$, and $x_1=a+h$, where $h=(b-a) / 2$. (See Figure 4.4)"

# ╔═╡ e7eed220-a52f-43d4-b431-57344875c0f7
let
	f1(x) = x^2
	f2(x) = x^4
	f3(x) = (x+1)^(-1)
	f4(x) = sqrt(1+x^2)
	f5(x) = sin(x)
	f6(x) = exp(x)
	xs = [0,2];
	xs2 = collect(0:0.1:2);
	
	part_a_trap = trapezoidal(0,2,f5,9)
	part_a_simp = simpson(xs,f5)
	part_a_simp2 = simpson(0,2,f5,4)
	Dict(
		:a_trap=>part_a_trap,
		:a_simp=>part_a_simp,
		:a_simp2=>part_a_simp2
	)
	# part_b_trap = trapezoidal(xs,x->f2.(x))
	# part_b_simp = simpson(xs,x->f2.(x))
	# Dict(
	# 	:a_trap=>part_b_trap,
	# 	:a_simp=>part_b_simp
	# )

end

# ╔═╡ 76715703-febf-454f-85b3-9bca7d8a44f4
md"## Measuring Precision"

# ╔═╡ 32bf866c-03ec-48d4-bf41-31586e248a92
md"# 4.4 Composite Numerical Integration"

# ╔═╡ 76bc523f-cbc4-49e0-9dbd-286b2a98a894
let
	I1 = simpson(0,4,exp)
	I2 = simpson(0,2,exp)+simpson(2,4,exp)
	I3 = sum(simpson(i,i+1,x->exp.(x)) for i in 0:3)
	I1,I2,I3,exp(4)-1
end

# ╔═╡ 4e814e85-56f6-4407-9de8-73a4c0d30eb6
let
	n1 = sqrt(π^3*1e5/(24))
	Iexact = cos(0)-cos(π)
	Itrap = trapezoidal(0,π,sin,250)
	abs(Iexact-Itrap), 2e-5
	n2 = (π^5*1e5/(360))^(1/4) 
	Isimp = simpson(0,π,sin,16)
	abs(Iexact-Isimp) ,0.00002
end

# ╔═╡ 0318bdb6-cec0-4f1f-bd0d-89c10553b211
md"# 5.1 The Elementary Theory of Initial-Value Problems"

# ╔═╡ a92e9f8f-06de-4176-aa1a-aed8a24bc71a
md"## Well-Posed Problems"

# ╔═╡ 6444c4a4-3015-46d1-8c22-8b28f54ad102
cm"""
- __Question__: How do we determine whether a particular problem has the property that small changes, or perturbations, in the statement of the problem introduce correspondingly small changes in the solution?
"""

# ╔═╡ 5189936d-6808-4d88-81da-454b992079eb
cm"""
The problem in $(eqref("five_three")) is called a __perturbed problem__ associated with the original proble.
"""

# ╔═╡ c6a776c3-c345-44d7-83f3-088f7eeb20b8
begin
	delta0_html  = @bind δ0 Slider(0:0.01:5, show_value=true)
	delta_html  = @bind δ Slider(0:0.01:9, show_value=true)
	cm"""
	<div style="display:flex;justify-content:space-between;">
	
	<div> 
	
	``\delta_0 =`` $(delta0_html) 
	
	</div>
	<!---
	<div> 
	
	``\delta=`` $delta_html 
	
	</div>
	--->
	</div>
	"""
end

# ╔═╡ 74d5b989-04e0-45db-a5c5-316541e04fd7
let
	f(t,y) = y-t^2 + 1 
	y0 = 0.5 + δ0
	problem1 = ODEProblem((y,p,t)->f(t,y),0.5,(0,2))
	problem2 = ODEProblem((y,p,t)->f(t,y),y0,(0,2))
	sol1 = solve(problem1)
	sol2 = solve(problem2)
	#Plot
	plot(sol1, linewidth = 2, title = "Problem: Example 3",
	    xaxis = "Time in seconds", yaxis = "y axis",
	    label = "Original Numerical Solution")
	plot!(sol2, linewidth = 2, 
	    xaxis = "Time in seconds", yaxis = "y axis",
	    label = "Numerical Solution Perturbed")
	plot!(sol1.t, t -> (t+1)^2 - 0.5exp(t), lw = 3, ls = :dash, label = "Analytical Solution")
	plot!(sol2.t, t -> (t+1)^2 - (0.5-δ0)exp(t), lw = 3, ls = :dash, label = "Analytical Solution For Perturbed")
	# plot!(sol.t, t -> (t+1)^2 + (δ0 +δ-0.5)exp(t) -δ, lw = 3, ls = :dash, label = "Analytical Solution for Perturbed")
end

# ╔═╡ 45df680f-d7a0-4112-9234-fdd69fa4f92c
cm"""
# 5.2 Euler's Method

- Euler's method is the most elementary approximation technique for solving initial-value problems. 
- Although it is seldom used in practice, the simplicity of its derivation can be used to illustrate the techniques involved in the construction of some of the more advanced techniques, without the cumbersome algebra that accompanies these constructions.
- The object of Euler's method is to obtain approximations to the well-posed initial-value problem
```math
\frac{d y}{d t}=f(t, y), \quad a \leq t \leq b, \quad y(a)=\alpha .
```

- A continuous approximation to the solution ``y(t)`` will not be obtained; instead, approximations to ``y`` will be generated at various values, called __mesh points__, in the interval ``[a, b]``. 
- Once the approximate solution is obtained at the points, the approximate solution at other points in the interval can be found by interpolation.
"""

# ╔═╡ 20ffd955-9b04-4562-b06b-8ee38a4abaa9
cm"""
Thus, __Euler's method__ is
```math
w_0=\alpha \text {, }
```
``w_{i+1}=w_i+h f\left(t_i, w_i\right), \quad`` for each ``i=0,1, \ldots, N-1``
"""

# ╔═╡ 842361ff-b144-463f-8ae1-9be8a2682723
function euler(f,t,α)
	h = t[2]-t[1]
	reduce((c,ti)->[c..., c[end]+h*f(ti,c[end])],t[1:end-1];init=[α])
end

# ╔═╡ a554a842-3b8d-45b9-9d40-466bbb879986
let
	f(t,y)= y - t^2 + 1
	N = 50
	h = 2/N
	tspan  = 0:h:2
	w = euler(f,tspan,0.5)
	plot(tspan,w, title="Example 1", label="Numerical Solution")
	plot!(t->(t+1)^2-0.5exp(t), label="Analytical Solution")
end

# ╔═╡ e41760c1-8b60-4efa-9bce-dd915fd9b671
let
	a,b = 0, 2
	L=1
	M = 0.5*exp(2)-2
	f(t,y)= y - t^2 + 1
	y(t) = (1+t)^2 -0.5exp(t)
	N = 10
	h = 2/N
	tspan  = 0:h:2
	w = euler(f,tspan,0.5)
	bound(ti) = (h*M/(2L)*(exp(L*(ti-a))-1))
	T = [
		"t_i" collect(tspan)'
		"actual" y.(tspan)'
		"approximation" w'
		"error" abs.(y.(tspan)'-w')
		"bound" reshape(bound.(tspan),1,N+1)
	]
	pretty_table(HTML,T,header=vcat("",map(i->"w$i",0:N)))
	
end

# ╔═╡ b77f106b-5943-4053-9afc-a91a1554781b
md"# 5.4 Runge-Kutta Methods"

# ╔═╡ bf206834-929c-43fd-a35d-9cd1aa2976b2
begin
	RK2(f,tspan,y0) = begin
		w = Vector{Float64}(undef,length(tspan))
		h = tspan[2]-tspan[1]
		h2 = h/2
		reduce((c,i)-> begin 
			wi = if i == 1
			 c 
			else 
				t = tspan[i-1]
				c+h*f(t+h2,c+h2*f(t,c))
			end
			w[i]=wi
			wi
		end,1:length(tspan);init=y0)
		w
	end
	MY_RK4(f,tspan,y0) = begin
		w = Vector{Float64}(undef,length(tspan))
		h = tspan[2]-tspan[1]
		h2 = h/2
		reduce((c,i)-> begin 
			wi = if i == 1
			 c 
			else 
				t = tspan[i-1]
				k1 = h*f(t,c)
				k2 = h*f(t+h2,c+0.5*k1)
				k3 = h*f(t+h2,c+0.5*k2)
				k4 = h*f(tspan[i],c+k3)
				
				twi = c + (1/6.0)*(k1+2k2+2k3+k4)
				
				twi
			end
			w[i]=wi
			wi
		end,1:length(tspan);init=y0)
		w
	end
	
end

# ╔═╡ 0c1b7f88-2907-474c-877b-4cfaf212dcf7
md"""## Runge-Kutta Methods of Order Two (Midpoint Method)
"""

# ╔═╡ a231aed6-741a-4fc7-b687-65042726dc3b
cm"""
```math
\begin{aligned}
w_0 & =\alpha \\
w_{i+1} & =w_i+h f\left(t_i+\frac{h}{2}, w_i+\frac{h}{2} f\left(t_i, w_i\right)\right), \quad \text { for } i=0,1, \ldots, N-1
\end{aligned}
```
"""

# ╔═╡ 99eefc06-0c8a-4326-af71-45aa94814703
let
	f(t,y) = y - t^2 + 1
	y0 = 0.5
	tspan = 0.0:0.2:2.0
	# w1 = RK2(f,tspan,y0)
	# problem = ODEProblem(f,y0,(0,2.0))
	w2 = MY_RK4(f,tspan,y0)
	# solve(problem,RK4())
end

# ╔═╡ 008d89d6-7cdc-4d37-bf36-50efed0d03be
md"## Higher-Order Runge-Kutta Methods"

# ╔═╡ c2169908-379d-4e5d-9f7f-d6fcecdf8f20
cm"""
### Runge-Kutta Order Four
```math
\begin{aligned}
w_0 & =\alpha, \\
k_1 & =h f\left(t_i, w_i\right), \\
k_2 & =h f\left(t_i+\frac{h}{2}, w_i+\frac{1}{2} k_1\right), \\
k_3 & =h f\left(t_i+\frac{h}{2}, w_i+\frac{1}{2} k_2\right), \\
k_4 & =h f\left(t_{i+1}, w_i+k_3\right), \\
w_{i+1} & =w_i+\frac{1}{6}\left(k_1+2 k_2+2 k_3+k_4\right),
\end{aligned}
```
"""

# ╔═╡ 4dd7bade-7523-4fa6-a862-25d2c61dbf9a
begin
	function add_space(n=1)
		repeat("&nbsp;",n)
	end
    function post_img(img::String, w=500)
        res = Resource(img, :width => w)
        cm"""
      <div class="img-container">

      $(res)

      </div>"""
    end
    function poolcode()
        cm"""
      <div class="img-container">

      $(Resource("https://www.dropbox.com/s/cat9ots4ausfzyc/qrcode_itempool.com_kfupm.png?raw=1",:width=>300))

      </div>"""
    end
    function define(t="")
        beginBlock("Definition", t)
    end
    function bbl(t)
        beginBlock(t, "")
    end
    function bbl(t, s)
        beginBlock(t, s)
    end
    ebl() = endBlock()
	function theorem(s)
		bth(s)
	end
    function bth(s)
        beginTheorem(s)
    end
    eth() = endTheorem()
    ex(n::Int; s::String="") = ex("Example $n", s)
    ex(t, s) = example(t, s)
    function beginBlock(title, subtitle)
        """<div style="box-sizing: border-box;">
       	<div style="display: flex;flex-direction: column;border: 6px solid rgba(200,200,200,0.5);box-sizing: border-box;">
       	<div style="display: flex;">
       	<div style="background-color: #FF9733;
       	    border-left: 10px solid #df7300;
       	    padding: 5px 10px;
       	    color: #fff!important;
       	    clear: left;
       	    margin-left: 0;font-size: 112%;
       	    line-height: 1.3;
       	    font-weight: 600;">$title</div>  <div style="olor: #000!important;
       	    margin: 0 0 20px 25px;
       	    float: none;
       	    clear: none;
       	    padding: 5px 0 0 0;
       	    margin: 0 0 0 20px;
       	    background-color: transparent;
       	    border: 0;
       	    overflow: hidden;
       	    min-width: 100px;font-weight: 600;
       	    line-height: 1.5;">$subtitle</div>
       	</div>
       	<p style="padding:5px;">
       """
    end
    function beginTheorem(subtitle)
        beginBlock("Theorem", subtitle)
    end
    function endBlock()
        """</p></div></div>"""
    end
    function endTheorem()
        endBlock()
    end
    function example(lable, desc)
        """<div style="display:flex;">
       <div style="
       font-size: 112%;
           line-height: 1.3;
           font-weight: 600;
           color: #f9ce4e;
           float: left;
           background-color: #5c5c5c;
           border-left: 10px solid #474546;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$lable:</div>
       <div style="flex-grow:3;
       line-height: 1.3;
           font-weight: 600;
           float: left;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$desc</div>
       </div>"""
    end
    @htl("")
end

# ╔═╡ 19d0bd5d-0168-4952-9ec3-3683424ce231
begin
    text_book = post_img("https://m.media-amazon.com/images/I/51ziKPbuEmL.jpg", 200)
    md""" # Syllabus
    ## Syallbus
    See here [Term 233 - MATH371 - Syllabus](https://www.dropbox.com/scl/fi/qxrcxxa1pxu3wctqzj0gg/T233_MATH371_Syllabus.pdf?rlkey=p715s0xldipiorxyfupe219og&raw=1)
    ## Textbook
    __Textbook: “Numerical Analysis” by Richard L. Burden, J. Douglas Faires 10th Edition (2016)__
    $text_book

    ## Office Hours
    I strongly encourage all students to make use of my office hours. These dedicated times are a valuable opportunity for you to ask questions, seek clarification on lecture material, discuss challenging problems, and get personalized feedback on your work. Engaging with me during office hours can greatly enhance your understanding of the course content and improve your performance. Whether you're struggling with a specific concept or simply want to delve deeper into the subject, I am here to support your learning journey. Don't hesitate to drop by; __your success is my priority__.

    | Day       | Time        |
    |-----------|-------------|
    | Monday    | 3:20-4:00PM |
    | Wednesday | 3:20-4:00PM |
    Also you can ask for an online meeting through __TEAMS__.
    """
end

# ╔═╡ b4c501ec-ed11-42ae-988e-6e73becf0d7e
cm"""

$(bth("Taylor Theorem"))

Suppose `` f \in C^n[a, b] ``, `` f^{(n+1)} `` exists on ``[a, b]``, and `` x_0 \in [a, b] ``. For every `` x \in [a, b] ``, there exists a number `` \xi(x) `` between `` x_0 `` and `` x `` with
```math
f(x) = P_n(x) + R_n(x),
```
where
```math
P_n(x) = f(x_0) + f'(x_0)(x - x_0) + \frac{f''(x_0)}{2!}(x - x_0)^2 + \cdots + \frac{f^{(n)}(x_0)}{n!}(x - x_0)^n
```
```math
= \sum_{k=0}^{n} \frac{f^{(k)}(x_0)}{k!}(x - x_0)^k
```
and

```math
R_n(x) = \frac{f^{(n+1)}(\xi(x))}{(n+1)!} (x - x_0)^{n+1}.
```
$(eth())
$(bbl("Remarks"))
- n``^{th}`` Taylor Polynomial ``P_n(x) ``
- Remainder Term `` R_n(x) ``
"""

# ╔═╡ 75115115-19c5-44c5-8c6a-7d3874228f35
cm"""
$(ex(1))
Let ``f(x) = \cos x`` and ``x_0 = 0``. Determine

- the second Taylor polynomial for ``f`` about ``x_0``; and
- the third Taylor polynomial for ``f`` about ``x_0``.

"""

# ╔═╡ c8beba83-e040-449a-b656-89b7daed7f7c
cm"""
$(bbl("Definition","Round-off Error"))
The error that is produced when a calculator or computer is used to perform real number calculations is called __round-off error__. 
$(ebl())

- It occurs because the arithmetic performed in a machine involves numbers with only a finite number of digits, with the result that calculations are performed with only approximate representations of the actual numbers. 
- In a computer, only a relatively small subset of the real number system is used for the representation of all the real numbers. This subset contains only rational numbers, both positive and negative, and stores the fractional part, together with an exponential part.
"""

# ╔═╡ 0feaa05f-9d53-48d1-b920-ea7c389105f2
cm"""
$(example("Example",""))
Consider the machine number

 	0 10000000011 10111001000100000000000000000000000000000000000000000

"""

# ╔═╡ edfa34e9-ee7d-4a0a-9025-a1857d8b9561
cm"""
$(ex(1)) 
Determine the five-digit 
- (a) chopping and 
- (b) rounding 

values of the irrational number ``\pi``.

"""

# ╔═╡ 244d084f-ed09-4d43-9377-8569995328a8
cm"""
$(bbl("Definition","") )
Suppose that ``p^*`` is an approximation to ``p``. 
- The actual error is ``p-p^*``, 
- the absolute error is ``\left|p-p^*\right|``, and 
- the relative error is 
```math
\frac{\left|p-p^*\right|}{|p|}, \text{ provided that } p \neq 0.
```
$(ebl())
"""

# ╔═╡ 321e0282-36ae-43dc-adcb-c7a35666a334
cm"""
$(ex(2)) 
Determine the actual, absolute, and relative errors when approximating ``p`` by ``p^*`` when
- (a) ``p=0.3000 \times 10^1`` and ``p^*=0.3100 \times 10^1``;
- (b) ``p=0.3000 \times 10^{-3}`` and ``p^*=0.3100 \times 10^{-3}``;
- (c) ``p=0.3000 \times 10^4`` and ``p^*=0.3100 \times 10^4``.
"""

# ╔═╡ e3f7eaa9-3fb2-4564-9e9e-961c2ce2e5ad
cm"""
 $(bbl("Definition","significant digits"))
 The number ``p^*`` is said to approximate ``p`` to ``t`` __significant digits__ (or figures) if ``t`` is the largest nonnegative integer for which
```math
\frac{\left|p-p^*\right|}{|p|} \leq 5 \times 10^{-t}
```
 $(ebl())
 """

# ╔═╡ f84cb603-95f5-42cc-afbc-c110e8e183c1
cm"""
$(example("Example",""))
```math
\begin{array}{l|l|l|l|l|l|l|l}p & 0.1 & 0.5 & 100 & 1000 & 5000 & 9990 & 10000\end{array}
```
What is the least upper bound for ``|p-p^*|`` if ``p^*`` agrees with ``p`` to __four__ significant figures.
"""

# ╔═╡ 892b86fe-9f31-4a38-9013-1f90e3b06389
cm"""
$(ex(3))
Suppose that ``x=\frac{5}{7}`` and ``y=\frac{1}{3}``. Use five-digit chopping for calculating ``x+y, x-y, x \times y``, and ``x \div y``.
"""

# ╔═╡ 966f78ae-65e1-4efa-b642-79f343e48a9b
cm"""
$(ex(4)) Suppose that in addition to ``x=\frac{5}{7}`` and ``y=\frac{1}{3}`` we have
```math
u=0.714251, \quad v=98765.9, \quad \text { and } \quad w=0.111111 \times 10^{-4},
```
so that
```math
f l(u)=0.71425 \times 10^0, \quad f l(v)=0.98765 \times 10^5, \quad \text { and } \quad f l(w)=0.11111 \times 10^{-4} .
```

Determine the five-digit chopping values of ``x \ominus u,(x \ominus u) \otimes w,(x \ominus u) \otimes v``, and ``u \oplus v``.
"""

# ╔═╡ c3720d44-e0a5-413a-aaf9-b11c6e6c442b
cm"""
$(example("Example",""))
Solve the following in 4-digit rounding arithmetic
```math
x^2 +62.10x+1=0
```
"""

# ╔═╡ e7fd2c69-097d-400f-aa20-19265ebdd2eb
cm"""
$(ex(6)) Evaluate ``f(x)=x^3-6.1 x^2+3.2 x+1.5`` at ``x=4.71`` using three-digit arithmetic.
"""

# ╔═╡ 3651b884-9583-49fb-95e8-269f349cf1ce
cm"""
$(define("Algorithm")) 
An __algorithm__ is a procedure that describes, in an unambiguous manner, a finite sequence of steps to be performed in a specified order. The object of the algorithm is to implement a procedure to solve a problem or approximate a solution to the problem.
$(ebl())

$(define("Pseudocode"))
A __pseudocode__ specifies the form of the input to be supplied and the form of the desired output. 
$(ebl())
- Not all numerical procedures give satisfactory output for arbitrarily chosen input. As a consequence, a stopping technique independent of the numerical technique is incorporated into each algorithm to avoid infinite loops.
- Two punctuation symbols are used in the algorithms:
  - A period (.) indicates the termination of a step.
  - A semicolon (;) separates tasks within a step.
"""

# ╔═╡ 053b3718-7c29-4a81-9cb8-c639e515aa07
cm"""
$(ex(1))
The ``N`` th Taylor polynomial for ``f(x)=\ln x`` expanded about ``x_0=1`` is
```math
P_N(x)=\sum_{i=1}^N \frac{(-1)^{i+1}}{i}(x-1)^i,
```
and the value of ``\ln 1.5`` to eight decimal places is 0.40546511 . Construct an algorithm to determine the minimal value of ``N`` required for
```math
\left|\ln 1.5-P_N(1.5)\right|<10^{-5}
```
without using the Taylor polynomial remainder term.
"""

# ╔═╡ f9e82af8-cbcb-4eab-8b66-227a5f34ccef
cm"""
$(bbl("Solution"))

__INPUT__ value ``x``, tolerance ``T O L``, maximum number of iterations ``M``. 

__OUTPUT__ degree ``N`` of the polynomial or a message of failure.

__Step 1__ Set ``N=1``;
```math
\begin{array}{ll}
& y=x-1 ; \\
& \text{ SUM }=0 ; \\
& \text{ POWER }=y ; \\
& \text{ TERM }=y ; \\
& \text{ SIGN }=-1 . \quad \text{(Used to implement alternation of signs.)}
\end{array}
```

__Step 2__ While ``N \leq M`` do Steps 3-5.

__Step 3__ Set SIGN ``=-`` SIGN; ``\quad`` (Alternate the signs.)
```math
\begin{aligned}
&\text{S U M}=\text{S U M}+\text{S I G N} \cdot \text{T E R M} ; \quad \text{(Accumulate the terms.)} \\
&\text{P O W E R}=\text{P O W E R} \cdot y ;\\
&\text{TERM} = \text{POWER}/(N+1). \quad \text{Calculate the next term.}
\end{aligned}
```

__Step 4__ If ``|\text{T E R M}|<\text{T O L}`` then (Test for accuracy.)

OUTPUT ( ``N`` );
STOP. (The procedure was successful.)

__Step 5__ Set ``N=N+1``. (Prepare for the next iteration. (End Step 2))

__Step 6__ OUTPUT ('Method Failed'); (The procedure was unsuccessful.) 

STOP.
"""

# ╔═╡ e8278f0e-74a6-4cc3-8164-269391162662
cm"""
$(define("1.17")) 
Suppose that ``E_0>0`` denotes an error introduced at some stage in the calculations and ``E_n`` represents the magnitude of the error after ``n`` subsequent operations.
- If ``E_n \approx C n E_0``, where ``C`` is a constant independent of ``n``, then the growth of error is said to be linear.
- If ``E_n \approx C^n E_0``, for some ``C>1``, then the growth of error is called exponential.
"""

# ╔═╡ a8848e7f-dfd2-438a-956e-daa91ac7ecad
cm"""
$(example("Example",""))
Consider the sequence ``\displaystyle p_n=\left(\frac{1}{3}\right)^n`` and use __five-digit rounding arithmatic__ to compute the terms of this sequence.
"""

# ╔═╡ a69d2160-3a9e-4038-9d87-5119e9e05466
cm"""
$(example("Example",""))
Now consider the sequence ``\displaystyle p_n=1-\frac{2}{3}n`` and use __five-digit rounding arithmatic__ to compute the terms of this sequence.
"""

# ╔═╡ 90bd960f-794d-496f-bf28-3b52abba90cc
cm"""
$(define("1.18")) 
Suppose ``\left\{\beta_n\right\}_{n=1}^{\infty}`` is a sequence known to converge to zero and ``\left\{\alpha_n\right\}_{n=1}^{\infty}`` converges to a number ``\alpha``. If a positive constant ``K`` exists with
```math
\left|\alpha_n-\alpha\right| \leq K\left|\beta_n\right|, \quad \text { for large } n,
```
then we say that ``\left\{\alpha_n\right\}_{n=1}^{\infty}`` converges to ``\alpha`` with rate, or order, of convergence ``O\left(\beta_n\right)``. (This expression is read "big oh of ``\beta_n`` ".) It is indicated by writing ``\alpha_n=\alpha+O\left(\beta_n\right)``.
"""

# ╔═╡ 62fb93e7-4acc-4466-b932-ce6cfdaf8d61
cm"""
$(ex(2))
Suppose that, for ``n \geq 1``,
```math
\alpha_n=\frac{n+1}{n^2} \quad \text { and } \quad \hat{\alpha}_n=\frac{n+3}{n^3} \text {. }
```
"""

# ╔═╡ 30833702-0444-4d4e-86c6-be97ae1456e7
cm"""
$(define("1.19"))
Suppose that ``\lim _{h \rightarrow 0} G(h)=0`` and ``\lim _{h \rightarrow 0} F(h)=L``. If a positive constant ``K`` exists with

```math
|F(h)-L| \leq K|G(h)|, \quad \text{for sufficiently small } h,
```
then we write ``F(h)=L+O(G(h))``.
"""

# ╔═╡ a3ad9375-b2bf-4b3b-9f84-d18d41265640
cm"""
$(ex(3))
Use the third Taylor polynomial about ``h=0`` to show that ``\cos h+\frac{1}{2} h^2=1+O\left(h^4\right)``.
"""


# ╔═╡ 9be1a640-c90f-4171-ab54-f6926dba25be
cm"""
$(ex(1))
Show that ``f(x)=x^3+4 x^2-10=0`` has a root in ``[1,2]`` and use the Bisection method to determine an approximation to the root that is accurate to at least within ``10^{-4}``.
"""

# ╔═╡ b11a67a2-b544-4996-9816-82bfbcd70a18
cm"""
$(bth("2.1"))
Suppose that ``f \in C[a, b]`` and ``f(a) \cdot f(b)<0``. The Bisection method generates a sequence ``\left\{p_n\right\}_{n=1}^{\infty}`` approximating a zero ``p`` of ``f`` with
```math
\left|p_n-p\right| \leq \frac{b-a}{2^n}, \quad \text { when } \quad n \geq 1
```
"""

# ╔═╡ e69ad96e-7796-4db4-9ae1-049ad0971f9c
cm"""
$(ex(2)) Determine the number of iterations necessary to solve ``f(x)=x^3+4 x^2-10=0`` with accuracy ``10^{-3}`` using ``a_1=1`` and ``b_1=2``.
"""

# ╔═╡ ddb60135-0438-4381-8284-053c464ec506
cm"""
$(define("2.2"))
The number ``p`` is a fixed point for a given function ``g`` if ``g(p)=p``.
"""

# ╔═╡ 668930ff-00d1-46b6-97a3-b27f5f6628c1
cm"""
$(ex(1)) Determine any fixed points of the function ``g(x)=x^2-2``.
"""

# ╔═╡ 9d2676de-fb2e-4a2c-8701-ee823bca5f71
cm"""
$(bth("2.3")) 
(i) If ``g \in C[a, b]`` and ``g(x) \in[a, b]`` for all ``x \in[a, b]``, then ``g`` has at least one fixed point in ``[a, b]``.
(ii) If, in addition, ``g^{\prime}(x)`` exists on ``(a, b)`` and a positive constant ``k<1`` exists with
```math
\left|g^{\prime}(x)\right| \leq k, \quad \text { for all } x \in(a, b),
```
then there is exactly one fixed point in ``[a, b]``. (See Figure 2.3.)
$(eth())
"""

# ╔═╡ 00594dfc-5a8b-4551-b300-bb52eee81e04
cm"""
$(post_img("https://www.dropbox.com/scl/fi/qdx5nwju090krjjywtw51/figure2.3.png?rlkey=sdlxcn8qvfywjyizmf2upctg2&raw=1"))
"""



# ╔═╡ a4e97910-11e5-4ad5-b9d2-7d5635e2990b
cm"""
$(ex(2)) 
Show that ``g(x)=\left(x^2-1\right) / 3`` has a unique fixed point on the interval ``[-1,1]``.
"""

# ╔═╡ b468696b-bb1b-44f6-8777-623b1b327b0b
cm"""
$(ex(3))
Consider 
```math
g(x)=3^{-x}\quad \text{ on } \quad [0,1].
```
Thoerem 2.3 does not guarantee the uniquness of the fixed point.
"""

# ╔═╡ a46e742e-9869-478c-b7a8-99267ceb9116
cm"""
$(post_img("https://www.dropbox.com/scl/fi/wwb8jccvt23artv0ky0j3/algorithm2.2_fixed_point.png?rlkey=lw78ogrp7skapnv2klsbox9pf&raw=1",700))
"""

# ╔═╡ b28b570c-44f3-49e9-9b94-eb6f2ed89bbf
cm"""
$(example("Example",""))
Solve ``x^3+4 x^2-10=0`` in the intervale ``[1,2]``.
"""

# ╔═╡ 12aa745e-9361-4af3-8c8b-7a2ffa83e874
cm"""
$(bth("2.4 (Fixed-Point Theorem)"))
Let ``g \in C[a, b]`` be such that ``g(x) \in[a, b]``, for all ``x`` in ``[a, b]``. Suppose, in addition, that ``g^{\prime}`` exists on ``(a, b)`` and that a constant ``0 < k <1 `` exists with
```math
\left|g^{\prime}(x)\right| \leq k, \quad \text { for all } x \in(a, b) .
```

Then, for any number ``p_0`` in ``[a, b]``, the sequence defined by
```math
p_n=g\left(p_{n-1}\right), \quad n \geq 1,
```
converges to the unique fixed point ``p`` in ``[a, b]``.
$(eth())
"""

# ╔═╡ e3a20e0f-1524-479f-83e8-6fc2093e320b
cm"""
$(bbl("Corollary", "2.5")) If ``g`` satisfies the hypotheses of Theorem 2.4 , then bounds for the error involved in using ``p_n`` to approximate ``p`` are given by
```math
\left|p_n-p\right| \leq k^n \max \left\{p_0-a, b-p_0\right\}
```
and
```math
\left|p_n-p\right| \leq \frac{k^n}{1-k}\left|p_1-p_0\right|, \quad \text { for all } \quad n \geq 1
```
$(ebl())
"""

# ╔═╡ 9291e2ac-56a3-41a0-88b8-d3d8ecb7e819
cm"""
- __Newton's method__ starts with an initial approximation __``p_0``__ and 
	- generates the sequence ``\left\{p_n\right\}_{n=0}^{\infty}``, by
```math
p_n=p_{n-1}-\frac{f\left(p_{n-1}\right)}{f^{\prime}\left(p_{n-1}\right)}, \quad \text { for } n \geq 1
```
See Figure 2.7

$(post_img("https://www.dropbox.com/scl/fi/ojynxcvtkpk8zta8whse7/fig2.7.png?rlkey=hk25cb8v4l7xv4705x9b3muu6&raw=1",700))
"""

# ╔═╡ 66b6ed61-d870-46d4-a5ad-cb6afec3a9dc
cm"""
$(post_img("https://www.dropbox.com/scl/fi/emcox2hdvabw08m0npeu0/algo2.3.png?rlkey=69o3qsf1tl2whxywgjd8eyvqc&raw=1",700))
"""

# ╔═╡ 2cbceb55-e953-4035-a0d3-e4ec236038ab
cm"""
$(ex(1))
Consider the function ``f(x)=\cos x-x=0``. Approximate a root of ``f`` using (a) a fixed-point method, and (b) Newton's method.
"""

# ╔═╡ e57d0aa5-661f-43e8-885c-fd0aa1ab4cc6
cm"""
$(bth("2.6")) Let ``f \in C^2[a, b]``. If ``p \in(a, b)`` such that ``f(p)=0`` and ``f^{\prime}(p) \neq 0``, then there exists a ``\delta>0`` such that Newton's method generates a sequence ``\left\{p_n\right\}_{n=1}^{\infty}`` converging to ``p`` for any initial approximation ``p_0 \in[p-\delta, p+\delta]``.
$(eth())
"""

# ╔═╡ 00e4aeb4-d28c-4de6-b298-a47a9d8ee3ab
cm"""
$(post_img("https://www.dropbox.com/scl/fi/c3f1pl27xv8ieo27lbbrg/fig2.9.png?rlkey=64rg0o9375jz0pum3ggd3amcm&raw=1",700))
"""

# ╔═╡ 4e765c64-eac2-4654-badf-222601c888b7
cm"""
$(post_img("https://www.dropbox.com/scl/fi/rw8czxi3l1gs2w50szffr/fig2.4.png?rlkey=0ruonucku445z3y4cp6xa84e8&raw=1",700))
"""

# ╔═╡ 079acdc9-781f-43c1-bd93-12e40142a0af
cm"""
$(ex(2))
Solve Example 1 using the __Secant Method__.
"""

# ╔═╡ fd29ef1c-d683-4a84-9595-04d0503d61ab
cm"""
- __Algebraic polynomials__, the set of functions of the form
```math
P_n(x)=a_n x^n+a_{n-1} x^{n-1}+\cdots+a_1 x+a_0,
```
``\color{white}{--.}``where ``n`` is a nonnegative integer and ``a_0, \ldots, a_n`` are real constants.
- Given __any function__, defined and continuous on a closed and bounded interval, __there exists a polynomial__ that is as "close" to the given function as desired. This result is expressed precisely in the Weierstrass Approximation Theorem. (See Figure 3.1.)

$(post_img("https://www.dropbox.com/scl/fi/uqu4r9frmxmrxa032yvae/fig3.1.png?rlkey=2l7ws8elwvptpkfu7omco1dqi&raw=1",700))


"""

# ╔═╡ 59929002-a335-4368-ab7c-a9ad080b3e78
cm"""
$(bth("3.1 (Weierstrass Approximation Theorem)")) 
Suppose ``f`` is defined and continuous on ``[a, b]``. For each ``\epsilon>0``, there exists a polynomial ``P(x)``, with the property that
```math
|f(x)-P(x)|<\epsilon, \quad \text { for all } x \text { in }[a, b] .
```
"""

# ╔═╡ eec05ea5-cb8a-4979-a5c9-23eebd9afe3f
cm"""
$(define("interpolation"))
Using a polynomial for approximation within the interval given by the endpoints is called __polynomial interpolation.__
"""

# ╔═╡ cc596f86-0762-4878-9c32-aaf3981b6398
cm"""
$(ex(1))
Determine the linear Lagrange interpolating polynomial that passes through the points ``(2,4)`` and ``(5,1)``.
"""

# ╔═╡ 54e544f0-f1d0-467d-9a7a-b773053895df
cm"""$(post_img("https://www.dropbox.com/scl/fi/tusc577dwl8j0vcgbubez/fig3.4.png?rlkey=osxczpcds3axenn0nzae7kxa0&raw=1",700))
"""

# ╔═╡ 1c83d413-5ab1-40c2-a533-25a312889f5c
cm"""$(post_img("https://www.dropbox.com/scl/fi/hgr57g6wf7np8hdn5un8x/fig3.5.png?rlkey=2v4sy729j19j33ok1lmnilukq&raw=1",700))"""

# ╔═╡ 1c4799ce-911b-4980-80fc-7c56c2b2a6ff
cm"""
$(ex(2))
- (a) Use the numbers (called nodes) ``x_0=2, x_1=2.75``, and ``x_2=4`` to find the second Lagrange interpolating polynomial for ``f(x)=1 / x``.
- (b) Use this polynomial to approximate ``f(3)=1 / 3``.
"""

# ╔═╡ 90297730-1473-4336-a884-d9441f3103a9
cm"""
$(bth("3.3"))
Suppose ``x_0, x_1, \ldots, x_n`` are distinct numbers in the interval ``[a, b]`` and ``f \in C^{n+1}[a, b]``. Then, for each ``x`` in ``[a, b]``, a number ``\xi(x)`` (generally unknown) between ``\min \left\{x_0, x_1, \ldots, x_n\right\}``, and the ``\max \left\{x_0, x_1, \ldots, x_n\right\}`` and hence in ``(a, b)``, exists with
```math
f(x)=P(x)+\frac{f^{(n+1)}(\xi(x))}{(n+1)!}\left(x-x_0\right)\left(x-x_1\right) \cdots\left(x-x_n\right),
```
where ``P(x)`` is the interpolating polynomial given above.
$(eth())
"""

# ╔═╡ e680937f-1fc5-4416-9991-8153bf604d64
cm"""
$(ex(3))
In Example 2, we found the second Lagrange polynomial for ``f(x)=1 / x`` on [2, 4] using the nodes ``x_0=2, x_1=2.75``, and ``x_2=4``. Determine the error form for this polynomial and the maximum error when the polynomial is used to approximate ``f(x)`` for ``x \in[2,4]``.
"""

# ╔═╡ ea5cb02b-49b4-4edf-89ab-d2fd6155e2a2
cm"""
$(bbl("",""))
So, ``P_n(x)`` can be rewritten in a form called __Newton's DividedDifference__:
```math
P_n(x)=f\left[x_0\right]+\sum_{k=1}^n f\left[x_0, x_1, \ldots, x_k\right]\left(x-x_0\right) \cdots\left(x-x_{k-1}\right)
```

The value of ``f\left[x_0, x_1, \ldots, x_k\right]`` is independent of the order of the numbers ``x_0, x_1, \ldots, x_k``,
$(ebl())
"""

# ╔═╡ c04f7b8a-f3ad-49d7-9610-48c5a7305649
cm"""
$(post_img("https://www.dropbox.com/scl/fi/36fy60w6547qrw1cu5mqu/algo3.2.png?rlkey=zyws01vns88ecijq3fww9cu5f&raw=1",700))
"""

# ╔═╡ 611bfd31-5eab-41ee-b410-120739748a2a
cm"""
$(ex(1)) Complete the divided difference table for the data used in the following Table and construct the interpolating polynomial that uses all these data.
| ``x`` | ``f(x)`` |
| :--- | :---: |
| 1.0 | 0.7651977 |
| 1.3 | 0.6200860 |
| 1.6 | 0.4554022 |
| 1.9 | 0.2818186 |
| 2.2 | 0.1103623 |
"""

# ╔═╡ 6a1d699c-fc16-472d-8395-a486890a089d
cm"""
The simplest piecewise-polynomial approximation is piecewise-linear interpolation, which consists of joining a set of data points
```math
\left\{\left(x_0, f\left(x_0\right)\right),\left(x_1, f\left(x_1\right)\right), \ldots,\left(x_n, f\left(x_n\right)\right)\right\}
```
by a series of straight lines, as shown in Figure 3.7.

$(post_img("https://www.dropbox.com/scl/fi/4sj9035slysb342sxqc5b/fig3.7.png?rlkey=2lh6s4qdqqt1euqtnc6mr03mq&raw=1",700))
"""

# ╔═╡ a9aa5b73-bde3-4907-971d-696e26ada195
cm"""$(post_img("https://www.dropbox.com/scl/fi/q9d37ya9w96xt43gg1vfi/fig3.8.png?rlkey=gw741eycxys78h6btem6pz72x&raw=1",700))"""

# ╔═╡ 47e27a13-4e60-40db-9ed7-5433e0230bd3
cm"""
$(define("3.10"))

Given a function ``f`` defined on ``[a, b]`` and a set of nodes ``a=x_0 < x_1 < \cdots < x_n=b``, a __cubic spline interpolant ``S`` for ``f``__ is a function that satisfies the following conditions:
- (a) ``S(x)`` is a cubic polynomial, denoted ``S_j(x)``, on the subinterval ``\left[x_j, x_{j+1}\right]`` for each ``j=0,1, \ldots, n-1``;
- (b) ``\quad S_j\left(x_j\right)=f\left(x_j\right)`` and ``S_j\left(x_{j+1}\right)=f\left(x_{j+1}\right)`` for each ``j=0,1, \ldots, n-1``;
- (c) ``S_{j+1}\left(x_{j+1}\right)=S_j\left(x_{j+1}\right)`` for each ``j=0,1, \ldots, n-2``; (Implied by (b).)
- (d) ``S_{j+1}^{\prime}\left(x_{j+1}\right)=S_j^{\prime}\left(x_{j+1}\right)`` for each ``j=0,1, \ldots, n-2``;
- (e) ``S_{j+1}^{\prime \prime}\left(x_{j+1}\right)=S_j^{\prime \prime}\left(x_{j+1}\right)`` for each ``j=0,1, \ldots, n-2``;
- (f) One of the following sets of boundary conditions is satisfied:
	- (i) ``S^{\prime \prime}\left(x_0\right)=S^{\prime \prime}\left(x_n\right)=0``
__(natural (or free) boundary)__;
	- (ii) ``S^{\prime}\left(x_0\right)=f^{\prime}\left(x_0\right)`` and ``S^{\prime}\left(x_n\right)=f^{\prime}\left(x_n\right) \quad`` __(clamped boundary)__.
	- (iii) ``S_0^{\prime\prime\prime}\left(x_1\right)=S_1^{\prime\prime\prime}\left(x_1\right)`` and ``S_{n-2}^{\prime\prime\prime}\left(x_{n-1}\right)=S_{n-1}^{\prime\prime\prime}\left(x_{n-1}\right)`` [or ``d_0=d_1`` and ``d_{n-2}=d_{n-1}``] __(Not-a-knot boundary)__.
	
"""

# ╔═╡ 1fd86b21-8457-4c1f-961c-81a105bda5a3
cm"""

$(bth("3.11"))
If ``f`` is defined at ``a=x_0 < x_1 < \cdots < x_n=b``, then ``f`` has a unique natural spline interpolant ``S`` on the nodes ``x_0, x_1, \ldots, x_n``; that is, a spline interpolant that satisfies the natural boundary conditions ``S^{\prime \prime}(a)=0`` and ``S^{\prime \prime}(b)=0``.
$(eth())
"""

# ╔═╡ cc599ac5-7000-45c9-bb5a-4529df717046
cm"""
$(post_img("https://www.dropbox.com/scl/fi/3d3jjo4a5uh2spth11dlm/algo3.4.png?rlkey=jvgva2sh5ysyif9d5y3sy1ez0&raw=1",700))
"""

# ╔═╡ 4177166a-d638-4739-91e7-d0c42a80392d
cm"""
$(ex(1)) Construct a natural cubic spline that passes through the points ``(1,2),(2,3)``, and ``(3,5)``.
"""


# ╔═╡ be8b89bb-e38f-423d-80b3-fa00fac7ad9a
cm"""
$(bbl("MATLAB",""))
We can use the MATLAB built-in function __`csape`__ (Cubic spline interpolation with end condition)


The `csape` function in MATLAB is used for cubic spline interpolation with specified end conditions. Here's a detailed explanation of the provided code snippet:


__Natural Boundary Conditions__
```matlab
xs = [1 2 3];
ys =[2 3 5];
pp = csape(xs,[0 ys 0],[2 2]);
```
##### Explanations:

1. **Vectors Definition**:
    - `xs = [1 2 3];` defines the `x` coordinates of the data points.
    - `ys =[2 3 5];` defines the corresponding `y` coordinates of the data points.

2. **csape Function Call**:
    - `pp = csape(xs,[0 ys' 0],[2 2]);` performs the cubic spline interpolation.
    - `csape` stands for "Cubic Spline with End Conditions".

##### Arguments:
- `xs`: The `x` coordinates of the data points.
- `[0 ys 0]`: The `y` coordinates are augmented by appending a `0` at the start and the end, `[0, ys', 0]`. 
    - The zeros (`0`) at the start and end are used to specify the boundary conditions.
- `[2 2]`: Specifies the boundary conditions type. 
    - The value `2` (second derivative) indicates that natural spline boundary conditions are used at both ends.
    - Natural spline boundary conditions ensure that the second derivatives at the end points are zero.

##### outputs
__`pp` Structure__ 
```{=matlab}
pp = 
  struct with fields:
    form: 'pp'
    breaks: [1 2 3]
    coefs: [2x4 double]
    pieces: 2
    order: 4
    dim: 1
```

The output `pp` is a piecewise polynomial structure that can be used for evaluating the spline at any desired point within the range of `xs`.

The `pp` structure is a MATLAB data type that contains information about a piecewise polynomial, including splines. It typically includes the following fields:

**form**: A string that specifies the form of the piecewise polynomial. For cubic splines, this is usually 'pp' (for piecewise polynomial). 
```{=matlab}
pp.form
``` 
This would output `'pp'`, indicating that this is a piecewise polynomial.

**breaks**: A vector of points where the pieces of the polynomial meet. These are the `x` coordinates of the data points provided for interpolation.
```{=matlab}
pp.breaks
```
This would output `[1 2 3]`, which are the x-coordinates of the original data points.

**coefs**: A matrix where each row contains the coefficients of the polynomial for a specific interval. For cubic splines, each row contains four coefficients corresponding to the cubic, quadratic, linear, and constant terms of the polynomial in that interval.
```{=matlab}
pp.coefs
```
This would output a 2x4 matrix, where each row contains the coefficients `[a, b, c, d]` of the polynomial for the corresponding interval. For example:
```{=matlab}
pp.coefs
ans =
	1.0000   -0.5000    0.5000    2.0000
	1.5000   -0.2500    0.7500    3.0000
```
Each row represents a cubic polynomial of the form:
```math
    p(x) = a(x - x_i)^3 + b(x - x_i)^2 + c(x - x_i) + d
```
for the corresponding interval ``[x_i, x_{i+1}]``.

**pieces**: The number of polynomial pieces, which is typically one less than the number of break points.
```{=matlab}
pp.pieces
```
This would output `2`, indicating there are 2 polynomial pieces.

**order**: The order of the polynomial. For cubic splines, this is 4 because the polynomial includes cubic, quadratic, linear, and constant terms.
```{=matlab}
pp.order
```
This would output `4`, indicating the polynomial is cubic (four coefficients).


**dim**: The dimension of the target. For univariate splines, this is 1.
```{=matlab}
pp.dim
```
This would output `1`, indicating the dimension of the target is univariate (one-dimensional).

##### Evaluating the Spline:

To evaluate the spline at a given point, you can use the `ppval` function:

```{=matlab}
x_eval = 1.5;
y_eval = ppval(pp, x_eval);
```


##### Example Usage:
To evaluate the spline at specific points, you can use the `ppval` function:

```{=matlab}
% Example evaluation points
evaluation_points = linspace(min(xs), max(xs), 100);

% Evaluate the spline at the desired points
spline_values = ppval(pp, evaluation_points);

% Plot the original data points and the interpolated spline
plot(xs, ys, 'o', evaluation_points, spline_values, '-');
legend('Data Points', 'Cubic Spline');
title('Cubic Spline Interpolation with Natural Boundary Conditions');
```

This script will plot the original data points along with the interpolated cubic spline, showing how the spline fits the data under natural boundary conditions.
$(ebl())
"""

# ╔═╡ 2482f975-0c63-45ab-8f9a-a2752809b192
cm"""
$(ex(2))
At the beginning of Chapter 3, we gave some Taylor polynomials to approximate the exponential ``f(x)=e^x``. Use the data points ``(0,1),(1, e),\left(2, e^2\right)``, and ``\left(3, e^3\right)`` to form a natural spline ``S(x)`` that approximates ``f(x)=e^x``.
"""

# ╔═╡ 755a1304-ec77-4ed8-806e-c287f13fdb89
cm"""
$(example("Example",""))
Approximate the integral of ``f(x)=e^x`` on ``[0,3]``
"""

# ╔═╡ 1e34cf47-8338-484e-b27a-fabfc7ef1b0b
cm"""
$(ex(3)) 
In Example 1, we found a natural spline ``S`` that passes through the points (1, 2), (2, 3), and ``(3,5)``. Construct a clamped spline ``s`` through these points that has ``s^{\prime}(1)=2`` and ``s^{\prime}(3)=1``.

"""

# ╔═╡ f0d5e130-d667-4f38-bfef-5aca795b880f
cm"""
$(bbl("MATLAB",""))
```{=matlab}
	xs =[1 2 3]
	ys = [2 3 5]
	% clamped with f'(1)=2 and f'(3)=1
 	pp = csape(xs,[2 ys 1],[1 1])
```
$(ebl())
"""

# ╔═╡ 766a993a-b41e-4de3-b511-6fa6e4b092e7
cm"""
$(ex(2))
Example 2 Fit the following data with the discrete least squares polynomial of degree at most 2 .
| ``i`` | ``x_i`` | ``y_i`` |
| :---: | :---: | :---: |
| 1 | 0 | 1.0000 |
| 2 | 0.25 | 1.2840 |
| 3 | 0.50 | 1.6487 |
| 4 | 0.75 | 2.1170 |
| 5 | 1.00 | 2.7183 |
"""

# ╔═╡ 3ade9674-04a8-4592-bb30-e807ffb0e17d
cm"""
$(ex(3)) Find exponential least squares for
| ``x_i`` | ``y_i`` |
| :---: | :---: |
| 1.00 | 5.10 |
| 1.25 | 5.79 |
| 1.50 | 6.53 |
| 1.75 | 7.45 |
| 2.00 | 8.46 |
"""

# ╔═╡ ac719200-b17f-4be9-be43-43b37b5c4018
cm"""
$(bbl("MATLAB",""))
The `fit` function in MATLAB is used to fit curves or surfaces to data. When using the `fit` function with different `fittype` options, you specify the type of model you want to fit to your data. Here is a brief explanation of using `fit` with the `fittype` options 'poly1', 'poly2', 'exp1', and 'power1':

##### Syntax
```{=matlab}
f = fit(xs, ys, 'fittype')
```

Here, `xs` and `ys` are vectors containing your data points, and `'fittype'` specifies the type of curve you want to fit.

##### `fittype` Options:

1. **'poly1'** (Linear Polynomial)
   - Fits a linear polynomial of the form `` f(x) = p1 \cdot x + p2 ``
   - **Example Usage**:
     ```{=matlab}
     f = fit(xs, ys, 'poly1');
     ```
   - **Description**: This fits a straight line to the data points, where `p1` is the slope and `p2` is the intercept.

2. **'poly2'** (Quadratic Polynomial)
   - Fits a quadratic polynomial of the form `` f(x) = p1 \cdot x^2 + p2 \cdot x + p3 ``
   - **Example Usage**:
     ```{=matlab}
     f = fit(xs, ys, 'poly2');
     ```
   - **Description**: This fits a parabola to the data points, where `p1`, `p2`, and `p3` are the coefficients of the quadratic, linear, and constant terms, respectively.

3. **'exp1'** (Single-Term Exponential)
   - Fits an exponential function of the form `` f(x) = a \cdot e^{b \cdot x} ``
   - **Example Usage**:
     ```{=matlab}
     f = fit(xs, ys, 'exp1');
     ```
   - **Description**: This fits an exponential curve to the data points, where `a` and `b` are the coefficients.

4. **'power1'** (Single-Term Power)
   - Fits a power function of the form `` f(x) = a \cdot x^b ``
   - **Example Usage**:
     ```{=matlab}
     f = fit(xs, ys, 'power1');
     ```
   - **Description**: This fits a power-law curve to the data points, where `a` is the coefficient and `b` is the exponent.

##### Example

Given data points in vectors `xs` and `ys`:
```{=matlab}
xs = [1; 2; 3; 4; 5];
ys = [2; 4; 6; 8; 10];
```

###### Linear Polynomial Fit:
```{=matlab}
f = fit(xs, ys, 'poly1');
```

###### Quadratic Polynomial Fit:
```{=matlab}
f = fit(xs, ys, 'poly2');
```

###### Exponential Fit:
```{=matlab}
f = fit(xs, ys, 'exp1');
```

###### Power Fit:
```{=matlab}
f = fit(xs, ys, 'power1');
```

In each case, the `fit` function will return a fit object `f` that contains the coefficients of the fitted model, which can be used to evaluate the model, plot it, or analyze its properties.
$(ebl())
"""

# ╔═╡ 34a1ae00-ae62-4782-bb03-9de7a8d26d1f
cm"""
$(ex(1)) Use the forward-difference formula to approximate the derivative of ``f(x)=\ln x`` at ``x_0=1.8`` using ``h=0.1, h=0.05``, and ``h=0.01`` and determine bounds for the approximation errors.
"""

# ╔═╡ bf7b6d08-4801-4757-87f0-c61e8966fec3
cm"""
```math
\displaystyle f^{\prime}\left(x_0\right)=\frac{1}{2 h}\left[-3 f\left(x_0\right)+4 f\left(x_0+h\right)-f\left(x_0+2 h\right)\right]+\frac{h^2}{3} f^{(3)}\left(\xi_0\right),
```

$(add_space(22)) where ``\xi_0`` lies between ``x_0`` and ``x_0+2 h``.

"""

# ╔═╡ 4f4fcb2c-a636-40e3-b82e-a2b8888ca7ae
cm"""
```math
f^{\prime}\left(x_0\right)=\frac{1}{2 h}\left[f\left(x_0+h\right)-f\left(x_0-h\right)\right]-\frac{h^2}{6} f^{(3)}\left(\xi_1\right),
```
$(add_space(22))where ``\xi_1`` lies between ``x_0-h`` and ``x_0+h``.
"""


# ╔═╡ b6163d7f-1447-461b-8912-f3c152fad7d5
cm"""
$(ex(2)) Values for ``f(x)=x e^x`` are given in below. Use all the applicable three-point formulas to approximate ``f^{\prime}(2.0)``.
| ``x`` | ``f(x)`` |
| :--- | :---: |
| 1.8 | 10.889365 |
| 1.9 | 12.703199 |
| 2.0 | 14.778112 |
| 2.1 | 17.148957 |
| 2.2 | 19.855030 |
"""

# ╔═╡ c6947376-99f5-4564-822c-0578f136b54a
cm"""
```math
f^{\prime \prime}\left(x_0\right)=\frac{1}{h^2}\left[f\left(x_0-h\right)-2 f\left(x_0\right)+f\left(x_0+h\right)\right]-\frac{h^2}{12} f^{(4)}(\xi),
```
$(add_space(15))for some ``\xi``, where ``x_0-h < \xi < x_0+h``.

$(add_space(10)) If ``f^{(4)}`` is continuous on ``\left[x_0-h, x_0+h\right]``, it is also bounded, and the approximation is ``O\left(h^2\right)``.
"""

# ╔═╡ f69d01bf-9832-4624-9c19-6c6eb024d24a
cm"""
$(ex(3)) Using the data from Example 2, use second derivative Midpoint formula to approximate ``f''(2.0)``
"""

# ╔═╡ c1894e87-0df3-431f-8ce6-f4dab38349c7
cm"""
$(bbl("MATLAB","<strong>diff</strong>"))

The `diff` function in MATLAB is used to calculate the difference between adjacent elements of an array. This function can be employed to approximate the first and second derivatives of a data set. Here is a brief explanation:

##### Syntax
```{=matlab}
Y = diff(X)
```

##### Description
- `Y = diff(X)` returns the differences between adjacent elements of `X`. If `X` is a vector, `Y` will be a vector of length `length(X)-1`. If `X` is a matrix, `Y` will be a matrix of the same number of rows, but with one less column, containing the differences between adjacent elements in each column.

##### Approximating the First Derivative
To approximate the first derivative of a function given by discrete data points, you can use the `diff` function to calculate the differences between consecutive data points.

### Example and Explanation

We will approximate the first and second derivatives of the sine function using finite differences.

#### Code:
```{=matlab}
h = 0.001;       % step size
X = -pi:h:pi;    % domain
f = sin(X);      % range
Y = diff(f)/h;   % first derivative
Z = diff(Y)/h;   % second derivative
plot(X(1:end-2), Z, 'k', X(1:end-1), Y, 'r', X, f, 'b')
legend('Second Derivative', 'First Derivative', 'Function')
```

#### Explanation:

1. **Step Size and Domain**:
    ```{=matlab}
    h = 0.001;       % step size
    X = -pi:h:pi;    % domain
    ```
    - `h` is the step size for the finite difference method.
    - `X` is a vector containing values from `-π` to `π` with increments of `h`.

2. **Function Values**:
    ```{=matlab}
    f = sin(X);      % range
    ```
    - `f` computes the sine of each value in `X`. This gives us the y-values of the sine function over the interval from `-π` to `π`.

3. **First Derivative Approximation**:
    ```{=matlab}
    Y = diff(f)/h;   % first derivative
    ```
    - `diff(f)` computes the difference between consecutive elements of `f`.
    - Dividing by `h` approximates the first derivative \( f'(x) \) using the central difference method.

4. **Second Derivative Approximation**:
    ```{=matlab}
    Z = diff(Y)/h;   % second derivative
    ```
    - `diff(Y)` computes the difference between consecutive elements of `Y`.
    - Dividing by `h` approximates the second derivative \( f''(x) \).

5. **Plotting the Results**:
    ```{=matlab}
    plot(X(1:end-2), Z, 'k', X(1:end-1), Y, 'r', X, f, 'b')
    legend('Second Derivative', 'First Derivative', 'Function')
    ```
    - `X(1:end-2)` is used for the second derivative because `Z` has two fewer elements than `X`.
    - `X(1:end-1)` is used for the first derivative because `Y` has one fewer element than `X`.
    - `X` is used for the original function `f`.
    - The plot shows the original sine function in blue, the first derivative (cosine) in red, and the second derivative (-sine) in black.

#### Explanation of Derivatives with the `diff` Function

##### First Derivative:
The `diff` function calculates the difference between consecutive elements in the vector `f`. This difference, divided by the step size `h`, gives an approximation of the first derivative using the finite difference method:

```math
Y_i = \frac{f_{i+1} - f_i}{h} 
```

This is essentially the slope of the line connecting two adjacent points on the sine curve, giving an estimate of the derivative at each point.

##### Second Derivative:
The second derivative is calculated by applying the `diff` function again to the first derivative `Y`. This provides the rate of change of the first derivative, which approximates the second derivative:

```math
Z_i = \frac{Y_{i+1} - Y_i}{h} = \frac{\left(\frac{f_{i+2} - f_{i+1}}{h}\right) - \left(\frac{f_{i+1} - f_i}{h}\right)}{h} 
```

Simplifying the above, we get:

```math
Z_i = \frac{f_{i+2} - 2f_{i+1} + f_i}{h^2} 
```

This represents the concavity or the rate of change of the slope, which is the second derivative of the function.

#### Conclusion

By using the `diff` function iteratively, we can approximate both the first and second derivatives of a function. The above example demonstrates this method for the sine function, showing how `diff` and finite differences can be used to estimate derivatives in MATLAB.

$(ebl())


"""

# ╔═╡ 979c7a6f-dd1d-4eb2-a862-7f370d654058
cm"""
- The need often arises for evaluating the definite integral of a function that has 
1. __no explicit antiderivative__ or 
2. whose __antiderivative is not easy to obtain__. 

- The basic method involved in approximating ``\int_a^b f(x) d x`` is called __numerical quadrature__. It uses a sum 
```math
\sum_{i=0}^n a_i f\left(x_i\right) \quad \text{to approximate}\quad \int_a^b f(x) d x.
```

- The methods of quadrature in this section are based on the interpolation polynomials given in Chapter 3. The basic idea is to select a set of distinct nodes ``\left\{x_0, \ldots, x_n\right\}`` from the interval ``[a, b]``. Then integrate the Lagrange interpolating polynomial
```math
P_n(x)=\sum_{i=0}^n f\left(x_i\right) L_i(x)
```
$(add_space(15))and its truncation error term over ``[a, b]`` to obtain
```math
\begin{aligned}
\int_a^b f(x) d x & =\int_a^b \sum_{i=0}^n f\left(x_i\right) L_i(x) d x+\int_a^b \prod_{i=0}^n\left(x-x_i\right) \frac{f^{(n+1)}(\xi(x))}{(n+1)!} d x \\
& =\sum_{i=0}^n a_i f\left(x_i\right)+\frac{1}{(n+1)!} \int_a^b \prod_{i=0}^n\left(x-x_i\right) f^{(n+1)}(\xi(x)) d x
\end{aligned}
```
"""

# ╔═╡ 6e52e8cf-5b78-41bd-b8fd-390e37614b76
cm"""
Consequently, 
```math
\begin{aligned}
\int_a^b f(x) d x & =\left[\frac{\left(x-x_1\right)^2}{2\left(x_0-x_1\right)} f\left(x_0\right)+\frac{\left(x-x_0\right)^2}{2\left(x_1-x_0\right)} f\left(x_1\right)\right]_{x_0}^{x_1}-\frac{h^3}{12} f^{\prime \prime}(\xi) \\
& =\frac{\left(x_1-x_0\right)}{2}\left[f\left(x_0\right)+f\left(x_1\right)\right]-\frac{h^3}{12} f^{\prime \prime}(\xi)
\end{aligned}
```
Using the notation ``h=x_1-x_0`` gives the following rule:
$(bbl("Trapezoidal Rule:",""))
```math
\int_a^b f(x) d x=\frac{h}{2}\left[f\left(x_0\right)+f\left(x_1\right)\right]-\frac{h^3}{12} f^{\prime \prime}(\xi) .
```
$(ebl())
"""

# ╔═╡ de6d1d86-94d3-4b8f-966c-4681dc4c5e2c
cm"""
$(post_img("https://www.dropbox.com/scl/fi/g1ejkkhsp9r0pya7buwzq/fig4.4.png?rlkey=zk35qvb4qfxxj0nqxib3ezjxy&raw=1",700))
"""

# ╔═╡ d426f020-40fe-4f10-a366-5ec99c11900b
cm"""
$(bbl("Simpson's Rule:",""))
```math
\int_{x_0}^{x_2} f(x) d x=\frac{h}{3}\left[f\left(x_0\right)+4 f\left(x_1\right)+f\left(x_2\right)\right]-\frac{h^5}{90} f^{(4)}(\xi) .
```
$(ebl())
"""

# ╔═╡ 00f5e3a3-dce4-4478-bf7c-2d1845205a71
cm"""
$(ex(1)) Compare the Trapezoidal rule and Simpson's rule approximations to ``\int_0^2 f(x) d x`` when ``f(x)`` is
- (a) ``x^2``
- (b) ``x^4``
- (c) ``(x+1)^{-1}``
- (d) ``\sqrt{1+x^2}``
- (e) ``\sin x``
- (f) ``e^x``
"""

# ╔═╡ dda247d8-e3b5-4d0f-a072-4a6f0dc262b5
cm"""
$(define("Degree of Accuracy"))
The degree of accuracy, or precision, of a quadrature formula is the largest positive integer ``n`` such that the formula is exact for ``x^k``, for each ``k=0,1, \ldots, n``.
"""

# ╔═╡ 8e77e17a-99a0-4c5d-8aa3-6e86f89f6292
cm"""
$(ex(1)) Use Simpson's rule to approximate ``\int_0^4 e^x d x`` and compare this to the results obtained by adding the Simpson's rule approximations for ``\int_0^2 e^x d x`` and ``\int_2^4 e^x d x`` and adding those for ``\int_0^1 e^x d x, \int_1^2 e^x d x, \int_2^3 e^x d x``, and ``\int_3^4 e^x d x``
"""

# ╔═╡ b086377d-2be1-485f-a241-12df11cde33e
cm"""
$(bth("4.4")) Let ``f \in C^4[a, b], n`` be even, ``h=(b-a) / n``, and ``x_j=a+j h``, for each ``j=0,1, \ldots, n``. There exists a ``\mu \in(a, b)`` for which the Composite Simpson's rule for ``n`` subintervals can be written with its error term as
```math
\int_a^b f(x) d x=\frac{h}{3}\left[f(a)+2 \sum_{j=1}^{(n / 2)-1} f\left(x_{2 j}\right)+4 \sum_{j=1}^{n / 2} f\left(x_{2 j-1}\right)+f(b)\right]-\frac{b-a}{180} h^4 f^{(4)}(\mu) \text {. }
```
$(eth())
"""

# ╔═╡ 4f1c3f5a-1490-4a1d-b26b-679773b90fac
cm"""
$(post_img("https://www.dropbox.com/scl/fi/ang29vx7rp1qly35kpss6/algo4.1.png?rlkey=g6b3yz0ew3habpwkvdut8vmld&raw=1",700))
"""

# ╔═╡ ab5d5e0d-afb3-457d-a169-e34608c1a93e
cm"""
$(bth("4.5"))
Let ``f \in C^2[a, b], h=(b-a) / n``, and ``x_j=a+j h``, for each ``j=0,1, \ldots, n``. There exists a ``\mu \in(a, b)`` for which the Composite Trapezoidal rule for ``n`` subintervals can be written with its error term as
```math
\int_a^b f(x) d x=\frac{h}{2}\left[f(a)+2 \sum_{j=1}^{n-1} f\left(x_j\right)+f(b)\right]-\frac{b-a}{12} h^2 f^{\prime \prime}(\mu)
```
"""

# ╔═╡ 17bcabbc-02c5-451e-9cbe-6e964cd306ef
cm"""
$(ex(2))
Determine values of ``h`` that will ensure an approximation error of less than 0.00002 when approximating ``\displaystyle\int_0^\pi \sin x d x`` and employing
<div style="display:flex;justify-content:space-evenly;">
<div>(a) Composite Trapezoidal rule and </div>
<div>(b) Composite Simpson's rule.</div>
</div>
"""

# ╔═╡ 07786494-0587-436d-a270-6f8e43d944e2
cm"""
$(bbl("MATLAB","integral"))
The MATLAB function `q = integral(fun,xmin,xmax)` is used to compute the numerical integral of a function over a specified interval. This function employs adaptive quadrature methods to evaluate the integral with high accuracy.

##### Syntax
```{=atlab}
q = integral(fun, xmin, xmax)
```

##### Description
- **`fun`**: A function handle representing the integrand, the function to be integrated. It should be of the form `fun = @(x) expression`.
- **`xmin`**: The lower limit of integration.
- **`xmax`**: The upper limit of integration.

##### Output
- **`q`**: The computed value of the integral over the interval `[xmin, xmax]`.

##### Example

Consider the following integral:
```math
\int_0^1 \sin(x) \, dx 
```

Let's compute this integral using `integral`.

###### Step-by-Step Solution

1. **Define the integrand function:**
    ```{=atlab}
    fun = @(x) sin(x);
    ```

2. **Set the limits of integration:**
    ```{=atlab}
    xmin = 0;
    xmax = 1;
    ```

3. **Call `integral` to compute the integral:**
    ```{=atlab}
    q = integral(fun, xmin, xmax);
    ```

4. **Display the result:**
    ```{=atlab}
    disp(['The integral of sin(x) from 0 to 1 is: ', num2str(q)])
    ```

##### Full Example Code
```{=atlab}
% Define the integrand function
fun = @(x) sin(x);

% Set the limits of integration
xmin = 0;
xmax = 1;

% Compute the integral
q = integral(fun, xmin, xmax);

% Display the result
disp(['The integral of sin(x) from 0 to 1 is: ', num2str(q)])
```

##### Explanation of the Example
- **`fun`**: This defines the integrand function \( \sin(x) \).
- **`xmin`**: The lower limit of integration is \( 0 \).
- **`xmax`**: The upper limit of integration is \( 1 \).
- **`integral`**: The function is called with the integrand, lower limit, and upper limit to compute the integral.
- **`disp`**: The result is displayed, showing the value of the integral.

##### Variations with Different Integrands

1. **Polynomial Integrand (fittype 'poly1'):**
    ```{=atlab}
    fun = @(x) x.^2;
    q = integral(fun, 0, 1);
    disp(['The integral of x^2 from 0 to 1 is: ', num2str(q)])
    ```

2. **Exponential Integrand (fittype 'exp1'):**
    ```{=atlab}
    fun = @(x) exp(x);
    q = integral(fun, 0, 1);
    disp(['The integral of exp(x) from 0 to 1 is: ', num2str(q)])
    ```

3. **Power Integrand (fittype 'power1'):**
    ```{=atlab}
    fun = @(x) x.^3;
    q = integral(fun, 0, 1);
    disp(['The integral of x^3 from 0 to 1 is: ', num2str(q)])
    ```

This example demonstrates how `integral` can be used to compute and display the value of an integral in MATLAB. By defining different integrand functions, you can use `integral` to compute various types of integrals over specified intervals.
$(ebl())
"""

# ╔═╡ 47d7efd0-48f6-479c-836a-f19c45d1318a
cm"""
$(define("5.1")) A function ``f(t, y)`` is said to satisfy a Lipschitz condition in the variable ``y`` on a set ``D \subset \mathbb{R}^2`` if a constant ``L>0`` exists with
```math
\left|f\left(t, y_1\right)-f\left(t, y_2,\right)\right| \leq L\left|y_1-y_2\right|,
```
whenever ``\left(t, y_1\right)`` and ``\left(t, y_2\right)`` are in ``D``. The constant ``L`` is called a Lipschitz constant for ``f``.
$(ebl())

$(ex(1)) Show that ``f(t, y)=t|y|`` satisfies a Lipschitz condition on the interval ``D=\{(t, y) \mid 1 \leq`` ``t \leq 2`` and ``-3 \leq y \leq 4\}``.
"""

# ╔═╡ ba9178c0-677b-470e-8164-7f54bd2438b5
cm"""
$(define("5.2"))
A set ``D \subset \mathbb{R}^2`` is said to be convex if whenever ``\left(t_1, y_1\right)`` and ``\left(t_2, y_2\right)`` belong to ``D``, then 
```math
\left((1-\lambda) t_1+\lambda t_2,(1-\lambda) y_1+\lambda y_2\right) \in D \quad\text{ for every }\lambda \in [0,1].
``` 
$(ebl())

$(post_img("https://www.dropbox.com/scl/fi/vtp12m6zn8yyz2tefl01s/fig5.1.png?rlkey=99r5x45z32zy51ozjfsqtkpdo&raw=1",700))
"""

# ╔═╡ 9ad7f65f-1e92-425c-9aa0-9804bbf1911d
cm"""
$(bth("5.3"))
Suppose ``f(t, y)`` is defined on a convex set ``D \subset \mathbb{R}^2``. If a constant ``L>0`` exists with
```math
\left|\frac{\partial f}{\partial y}(t, y)\right| \leq L, \quad \text { for all }(t, y) \in D,
```
then ``f`` satisfies a Lipschitz condition on ``D`` in the variable ``y`` with Lipschitz constant ``L``.
$(eth())
"""

# ╔═╡ e60aaf53-b008-40aa-8403-1f619b81e5b7
cm"""
$(theorem("5.4")) Suppose that ``D=\{(t, y) \mid a \leq t \leq b`` and ``-\infty < y < \infty\}`` and that ``f(t, y)`` is continuous on ``D``. If ``f`` satisfies a Lipschitz condition on ``D`` in the variable ``y``, then the initial-value problem
```math
y^{\prime}(t)=f(t, y), \quad a \leq t \leq b, \quad y(a)=\alpha,
```
has a unique solution ``y(t)`` for ``a \leq t \leq b``.
$(eth())

$(ex(2)) Use Theorem 5.4 to show that there is a unique solution to the initial-value problem
```math
y^{\prime}=1+t \sin (t y), \quad 0 \leq t \leq 2, \quad y(0)=0
```
"""

# ╔═╡ 427b3198-6485-462b-8d13-23322c96a8ac
cm"""
$(bbl("Exercise",""))
Use Theorem 5.4 to show that  the following initial-value problems has a unique solution and find the solution.
```math
y^{\prime}=-\frac{2}{t} y+t^2 e^t, 1 \leq t \leq 2, y(1)=\sqrt{2} e.
```
"""

# ╔═╡ 0a5b7cc2-a84e-403e-bc7d-b40f94550383
cm"""
$(define("5.5")) The initial-value problem
```math
\frac{d y}{d t}=f(t, y), \quad a \leq t \leq b, \quad y(a)=\alpha,
```
is said to be a __well-posed problem__ if:
- A unique solution, ``y(t)``, to the problem exists, and
- There exist constants ``\varepsilon_0>0`` and ``k>0`` such that for any ``\varepsilon``, in ``\left(0, \varepsilon_0\right)``, whenever ``\delta(t)`` is continuous with ``|\delta(t)|<\varepsilon`` for all ``t`` in ``[a, b]``, and when ``\left|\delta_0\right|<\varepsilon``, the initial-value problem
$(texeq"
\frac{d z}{d t}=f(t, z)+\delta(t), \quad a \leq t \leq b, \quad z(a)=\alpha+\delta_0,
\label{five_three}")
$(add_space(10))has a unique solution ``z(t)`` that satisfies
```math
|z(t)-y(t)| < k \varepsilon \quad \text { for all } t \text { in }[a, b]
```
"""

# ╔═╡ 4c56c9df-9186-4369-9021-c049a2eefdd0
cm"""
$(theorem("5.6"))
Suppose ``D=\{(t, y) \mid a \leq t \leq b`` and ``-\infty < y <\infty\}``. If ``f`` is continuous and satisfies a Lipschitz condition in the variable ``y`` on the set ``D``, then the initial-value problem
```math
\frac{d y}{d t}=f(t, y), \quad a \leq t \leq b, \quad y(a)=\alpha
```
is well posed.
$(eth())

$(ex(3)) Show that the initial-value problem
```math
\frac{d y}{d t}=y-t^2+1, \quad 0 \leq t \leq 2, \quad y(0)=0.5,
```
is well posed on ``D=\{(t, y) \mid 0 \leq t \leq 2`` and ``-\infty < y <\infty\}``.
"""

# ╔═╡ 3d590d35-f866-404d-9c2a-e8dc17712113
cm"""
$(bbl("Illustrations",""))
As an illustration, consider the solution to the perturbed problem
```math
\frac{d z}{d t}=z-t^2+1+\delta, \quad 0 \leq t \leq 2, \quad z(0)=0.5+\delta_0
```
$(ebl())
"""

# ╔═╡ cf37eeeb-7317-4050-8218-2a67a2432852
cm"""
$(post_img("https://www.dropbox.com/scl/fi/lug6uovsoetnx45f5dwwt/algo5.1.png?rlkey=p2r5n84k6jhp4nts5mvrah03h&raw=1",700))
"""

# ╔═╡ 8566080c-a6bc-4823-b1f2-747687cdada1
cm"""
$(ex(1)) Euler's method was used in the first illustration with ``h=0.5`` to approximate the solution to the initial-value problem
```math
y^{\prime}=y-t^2+1, \quad 0 \leq t \leq 2, \quad y(0)=0.5 .
```

Use Algorithm 5.1 with ``N=10`` to determine approximations and compare these with the exact values given by ``y(t)=(t+1)^2-0.5 e^t``.
"""

# ╔═╡ 66820a81-3ad6-4905-9810-0ab5dd38ce69
cm"""
$(theorem("5.9")) Suppose ``f`` is continuous and satisfies a Lipschitz condition with constant ``L`` on
```math
D=\{(t, y) \mid a \leq t \leq b \text { and }-\infty < y< \infty\}
```
and that a constant ``M`` exists with
```math
\left|y^{\prime \prime}(t)\right| \leq M, \quad \text { for all } t \in[a, b]
```
where ``y(t)`` denotes the unique solution to the initial-value problem
```math
y^{\prime}=f(t, y), \quad a \leq t \leq b, \quad y(a)=\alpha .
```

Let ``w_0, w_1, \ldots, w_N`` be the approximations generated by Euler's method for some positive integer ``N``. Then, for each ``i=0,1,2, \ldots, N``,
```math
\left|y\left(t_i\right)-w_i\right| \leq \frac{h M}{2 L}\left[e^{L\left(t_i-a\right)}-1\right]
```
"""

# ╔═╡ 91cf0689-7e38-4c5d-b718-b444c73853f6
cm"""
$(ex(2)) The solution to the initial-value problem
```math
y^{\prime}=y-t^2+1, \quad 0 \leq t \leq 2, \quad y(0)=0.5,
```
was approximated in Example 1 using Euler's method with ``h=0.2``. Use the inequality in Theorem 5.9 to find a bound for the approximation errors and compare these to the actual errors.
"""

# ╔═╡ fc5b69cd-eddd-4f1a-92c4-df4d26939821
cm"""
$(ex(2)) Use the Midpoint method  with ``N=10, h=0.2, t_i=0.2 i``, and ``w_0=0.5`` to approximate the solution to our usual example,
```math
y^{\prime}=y-t^2+1, \quad 0 \leq t \leq 2, \quad y(0)=0.5
```
"""

# ╔═╡ 767b00e9-3d6c-413c-8a82-1e9152aa5052
cm"""
$(bbl("MATLAB","ode45"))
The MATLAB function `[t,y] = ode45(odefun,tspan,y0)` is used to solve ordinary differential equations (ODEs) numerically using the Dormand-Prince method, which is a method of the Runge-Kutta family. This function is particularly useful for solving non-stiff differential equations.

##### Syntax
```{=matlab}
[t,y] = ode45(odefun, tspan, y0)
```

##### Description
- **`odefun`**: A function handle that defines the differential equation to be solved. It should be of the form `dy/dt = f(t,y)`.
- **`tspan`**: A vector specifying the interval of integration `[t0 tf]`. The solver will return the solution at each time step within this interval.
- **`y0`**: A vector specifying the initial conditions for the ODE.

##### Outputs
- **`t`**: A column vector of time points at which the solution was evaluated.
- **`y`**: A matrix where each row corresponds to the solution at a time point in `t`.

##### Example

Consider the following differential equation:
```math
\frac{dy}{dt} = -2y + t 
```

Let's solve this ODE using `ode45` over the interval [0, 5] with the initial condition `` y(0) = 1 ``.

###### Step-by-Step Solution

1. **Define the ODE function:**
    ```{=matlab}
    function dydt = odefun(t, y)
        dydt = -2*y + t;
    end
    ```

2. **Set the time span and initial condition:**
    ```{=matlab}
    tspan = [0 5];
    y0 = 1;
    ```

3. **Call `ode45` to solve the ODE:**
    ```{=matlab}
    [t, y] = ode45(@odefun, tspan, y0);
    ```

4. **Plot the results:**
    ```{=matlab}
    plot(t, y)
    xlabel('Time t')
    ylabel('Solution y')
    title('Solution of dy/dt = -2y + t using ode45')
    ```

##### Full Example Code
```{=matlab}
% Define the ODE function
function dydt = odefun(t, y)
    dydt = -2*y + t;
end

% Set the time span and initial condition
tspan = [0 5];
y0 = 1;

% Solve the ODE
[t, y] = ode45(@odefun, tspan, y0);

% Plot the results
plot(t, y)
xlabel('Time t')
ylabel('Solution y')
title('Solution of dy/dt = -2y + t using ode45')
```

##### Explanation of the Example
- **`odefun`**: This function defines the ODE `` \frac{dy}{dt} = -2y + t ``.
- **`tspan`**: The solver integrates the ODE from `` t = 0 `` to `` t = 5 ``.
- **`y0`**: The initial condition is `` y(0) = 1 ``.
- **`ode45`**: The function is called with the ODE function handle, time span, and initial condition to compute the solution.
- **`plot`**: The solution is plotted with time on the x-axis and the solution on the y-axis.

This example demonstrates how `ode45` can be used to solve and visualize the solution of an ODE in MATLAB.


$(ebl())

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Integrals = "de52edbc-65ea-441a-8357-d3a637375a31"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NonlinearSolve = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QRCoders = "f42e9828-16f3-11ed-2883-9126170b272d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "bee9ac6c13ea91562613eafd217676abd3afe64e"

[[deps.ADTypes]]
git-tree-sha1 = "2314e58e823f0fd6ee02dbbecb997370f501dd4a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "c0d491ef0b135fd7d63cbc6404286bc633329425"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.36"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ed2ec3c9b483842ae59cd273834e5b46206d6dda"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.11.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "420e2853770f50e5306b9d96b5a66f26e7c73bc6"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.9.4"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "585a387a490f1c4bd88be67eea15b93da5e85db7"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "2c6b7bf16fd850c551a765e313e7522ba455cbfd"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.151.4"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "Compat", "DocStringExtensions", "FillArrays", "LinearAlgebra", "PackageExtensionCompat", "SparseArrays", "SparseMatrixColorings"]
git-tree-sha1 = "23c6df13ad8fcffde4b0596d798911d2e309fc2c"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.5.5"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTapirExt = "Tapir"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tapir = "07d77754-e150-4737-8c94-cd238a1fb45b"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "490392af2c7d63183bfa2c8aaa6ab981c5ba7561"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.14"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "30a1848c4f4fc35d1d4bbbd125650f6a11b5bc6c"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.7"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "3a3177ba05b4763234819060fb6c2e1613379ca6"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.7.6"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "94997910aca72897524d2237c41eb852153b0f65"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "2be93e36303143c6fffd07e2222bbade35194d9e"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.3"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2de436b72c3422940cbe1367611d137008af7ec3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "3e527447a45901ea392fe12120783ad6ec222803"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.6"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "182c478a179b267dd7a741b6f8f4c3e0803795d6"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.6+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "801aef8228f7f04972e596b09d4dba481807c913"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.4"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "334d300809ae0a68ceee3444c6e99ded412bf0b3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.11.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "10f37537bbd83e52c63abf6393f209dbd641fedc"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.6.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d65554bad8b16d9562050c67e7223abf91eaba2f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.13+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.Integrals]]
deps = ["CommonSolve", "HCubature", "LinearAlgebra", "MonteCarloIntegration", "QuadGK", "Reexport", "SciMLBase"]
git-tree-sha1 = "ebf5737d823873add85809f2b52e20e3eae71997"
uuid = "de52edbc-65ea-441a-8357-d3a637375a31"
version = "4.4.1"

    [deps.Integrals.extensions]
    IntegralsArblibExt = "Arblib"
    IntegralsCubaExt = "Cuba"
    IntegralsCubatureExt = "Cubature"
    IntegralsFastGaussQuadratureExt = "FastGaussQuadrature"
    IntegralsForwardDiffExt = "ForwardDiff"
    IntegralsMCIntegrationExt = "MCIntegration"
    IntegralsZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.Integrals.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Cuba = "8a292aeb-7a57-582c-b821-06e4c11590b1"
    Cubature = "667455a9-e2ce-5579-9412-b964f529a492"
    FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MCIntegration = "ea1e2de9-7db7-4b42-91ee-0cd1bf6df167"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be50fe8df3acbffa0274a744f1a99d29c45a57f4"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "e7cbed5032c4c397a6ac23d1493f3289e01231c4"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.14"
weakdeps = ["Dates"]

    [deps.InverseFunctions.extensions]
    DatesExt = "Dates"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "267dad6b4b7b5d529c76d40ff48d33f7e94cb834"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "e459fda6b68ea8684b3fcd513d2fd1e5130c4402"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.16.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "fa518725cf83bf6f0c1aecd06adb2ea2d94c4dc6"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.0.5"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "7648cc20100504f4b453917aacc8520e9c0ecfb3"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.30.1"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = ["Enzyme", "EnzymeCore"]
    LinearSolveFastAlmostBandedMatricesExt = ["FastAlmostBandedMatrices"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8f6786d8b2b3248d79db3ad359ce95382d5a6df8"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.170"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "27d162f37cc29de047b527dab11a826dd3a650ad"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "1b9e613f2ca3b6cdcbfe36381e17ca2b66d4b3a1"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.3"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MonteCarloIntegration]]
deps = ["Distributions", "QuasiMonteCarlo", "Random"]
git-tree-sha1 = "722ad522068d31954b4a976b66a26aeccbf509ed"
uuid = "4886b29c-78c9-11e9-0a6e-41e1f4161f7b"
version = "0.2.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "5c1d1d9361e1417e5a065e1f84dc3686cbdaea21"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "898c56fbf8bf71afb0c02146ef26f3a454e88873"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.5"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "40325dcea1cb84a108efe05966bbb1f4b98e5eea"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.13.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSymbolicsExt = "Symbolics"
    NonlinearSolveZygoteExt = "Zygote"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "f4cb457ffac5f5cf695699f82c537073958a6a6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.2+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b4cde20f0e8c67fd35863794d5e548722f7bb71d"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.84.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoDevMacros]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "c3839362a712e6d9c2845d179edafe74371cb77b"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.7.4"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "93d8c75734da9192d0639406fe6fb446be0fba4f"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "b3e2bae88cf07baf0a051fe09666b8ef97aefe93"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.14"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "406c29a7f46706d379a3bce45671b4e3a39ddfbc"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.22"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "763a8ceb07833dd51bb9e3bbca372de32c0605ad"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.0"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QRCoders]]
deps = ["FileIO", "ImageCore", "ImageIO", "ImageMagick", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "b3e5fcc7a7ade2d43f0ffd178c299b7a264c268a"
uuid = "f42e9828-16f3-11ed-2883-9126170b272d"
version = "1.4.5"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.QuasiMonteCarlo]]
deps = ["Accessors", "ConcreteStructs", "LatticeRules", "LinearAlgebra", "Primes", "Random", "Requires", "Sobol", "StatsBase"]
git-tree-sha1 = "cc086f8485bce77b6187141e1413c3b55f9a4341"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.3.3"
weakdeps = ["Distributions"]

    [deps.QuasiMonteCarlo.extensions]
    QuasiMonteCarloDistributionsExt = "Distributions"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "3400ce27995422fb88ffcd3af9945565aad947f0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.23.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d483cd324ce5cf5d61b77930f0bbd6cb61927d21"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.2+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "7a6c5c8c38d2e37f45d4686c3598c20c1aebf48e"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.41.3"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools", "Setfield", "SparseArrays", "StaticArraysCore"]
git-tree-sha1 = "10499f619ef6e890f3f4a38914481cc868689cd5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.8"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "6ab4beaf88dcdd2639bead916f2347f81dcacd0e"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "913754ccbbc78720a4542b56a6bdfbab1c84c8f2"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.10.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "469f51f8c4741ce944be2c0b65423b518b1405b0"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.19.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "Compat", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "eed2446b3c3dd58f6ded3168998b8b2cb3fc9229"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.3.3"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "6e00379a24597be4ae1ee6b2d882e15392040132"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.5"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "25349bf8f63aa36acbff5e3550a86e9f5b0ef682"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.6"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "a5f6f138b740c9d93d76f0feddd3092e6ef002b7"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.22"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fb099adbd7504f1e68b4512828e9d94197a8b889"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "79813de27af70906d223fbd89ad90dba3d88a2b0"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "2.0.2"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "ForwardDiff", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "616ba6333a9a8132108ad2c283c53682860e58e9"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.30.4"

    [deps.Symbolics.extensions]
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxCoreExt = "LuxCore"
    SymbolicsPreallocationToolsExt = "PreallocationTools"
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LuxCore = "bb33d45b-7691-41d6-9220-0943567d0623"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "6f0cee95e74d1f6891ba6b35b8b219fd3d11b567"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.4.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[deps.TranscodingStreams]]
git-tree-sha1 = "d73336d81cafdc277ff45558bb7eaa2b04a8e472"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.10"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "66c68a20907800c0b7c04ff8a6164115e8747de2"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.0"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeTypeAbstraction", "LazyModules", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "ae67ab0505b9453655f7d5ea65183a1cd1b3cfa0"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.12.4"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e863582a41c5731f51fd050563ae91eb33cf09be"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.68"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─8ca0d1c5-166d-44f0-a17e-a6207c19459a
# ╟─dbc554cb-3061-48b8-aaa2-c6c13e7a3d75
# ╟─21331300-b8cb-410f-88a3-1feb63e55659
# ╟─a279f75a-1d5a-4c0e-aae6-a168dda0a277
# ╟─f581c560-afc6-460d-bb7e-48e13c1909d4
# ╟─7d454fff-b638-4678-92f6-12d816b541b1
# ╟─830e0fb0-8bfa-48b6-9210-cc0dc009d042
# ╟─19d0bd5d-0168-4952-9ec3-3683424ce231
# ╟─d7267eec-8b39-454e-b1ec-c7ae580190c4
# ╟─b4c501ec-ed11-42ae-988e-6e73becf0d7e
# ╠═cb26e993-d7c1-4e69-82f3-dcb20d1a4f37
# ╟─75115115-19c5-44c5-8c6a-7d3874228f35
# ╠═b54035ab-2813-4cdc-887a-e11625ede4aa
# ╠═e223e89c-af0c-44ed-bb83-fd51124c7899
# ╠═9e344be4-e96a-4b1f-9bac-52deb66c2658
# ╟─d60047cf-47fd-406f-b100-be2b1d195612
# ╟─7e6509af-4072-45f4-a3cf-d2891c586d58
# ╟─ee13d4c2-edce-4d7c-8ace-103aab1f7ba0
# ╟─df80d1bc-1285-48f3-8de1-5bf9631f078c
# ╟─afbd4725-9c76-4453-a62a-ebac185655a2
# ╟─422a3ccc-1e62-407b-bd3e-0703f1a5a33d
# ╟─2d460985-77a6-4453-b639-3312c98b742f
# ╟─9b0cd30c-3162-4423-b32f-88c9a01c1cbd
# ╟─a4214d60-b66e-4819-84ad-37653e45184e
# ╟─f610361d-e308-49be-8100-ad8b3882bac0
# ╠═98652338-3d5c-40f9-9220-c4f6e94f532c
# ╟─c8beba83-e040-449a-b656-89b7daed7f7c
# ╟─e4a92b5e-18a7-4007-aab3-3eeca447b8de
# ╠═726a42a1-9746-4083-88f3-0661cd7396b5
# ╟─11c16f32-f86b-43ba-8a58-532609c226e0
# ╟─0feaa05f-9d53-48d1-b920-ea7c389105f2
# ╠═4d8df2cc-74b1-445e-9cab-cb58af4d56bb
# ╠═2563dd79-3e17-4b81-82f3-43d6aea85b71
# ╠═c9819b7f-4bf8-46e3-a4bc-e991bf26ce7a
# ╠═13ac7d12-0740-45df-831c-30670afd9f40
# ╠═ed835d55-deaf-44b5-b360-f8eff7800f6d
# ╠═239e7eaa-b8da-47d4-b642-8415ef1e0683
# ╠═36211411-8f83-4373-8dd0-7b7f334f8bc1
# ╠═bb4ed82e-3b7e-4898-987f-9a9bac49c205
# ╠═66fa9c4b-f2f3-44a1-8755-bc13ffccaba1
# ╟─75bac38d-3c93-4cad-8302-6eb60e40038b
# ╟─ff2fe07a-d35b-4d5b-bf95-f669fafc2132
# ╟─9ddc2701-bf18-4f96-9394-36a95b399252
# ╟─edfa34e9-ee7d-4a0a-9025-a1857d8b9561
# ╠═7caaa18e-8ef8-4ae5-bd1c-e9fd0536b901
# ╟─43e1134e-020c-4fff-8164-050f0fcb9c38
# ╟─244d084f-ed09-4d43-9377-8569995328a8
# ╟─321e0282-36ae-43dc-adcb-c7a35666a334
# ╟─8c6005c5-e947-4a8a-8083-78ab5f84f07b
# ╠═c8d040f0-3264-42c6-84c5-be8c633c36e4
# ╟─e3f7eaa9-3fb2-4564-9e9e-961c2ce2e5ad
# ╟─f84cb603-95f5-42cc-afbc-c110e8e183c1
# ╟─d38e3f95-71c4-45dc-8ab9-0a216de91311
# ╟─786351e5-47cc-4d61-a47d-771bc47aa3a9
# ╟─2e3eefec-a29b-4dff-a863-eabbb4a3c7ef
# ╟─712e7cbd-1280-4885-8adb-d83de2c26560
# ╟─4db7b813-d0da-4fc6-8447-aae94fbae670
# ╟─892b86fe-9f31-4a38-9013-1f90e3b06389
# ╟─9010d6bd-26d2-4211-ab12-ac58d77edbbc
# ╟─966f78ae-65e1-4efa-b642-79f343e48a9b
# ╟─b50497e4-e0f9-4da6-bba3-48628cee0173
# ╟─c3720d44-e0a5-413a-aaf9-b11c6e6c442b
# ╠═379608bf-906d-45db-9cab-018a50fcbd07
# ╟─1922d625-ba5e-4ea2-b050-a3c64077c742
# ╟─5ed7d51c-b5fa-4953-b2b4-c2040edced33
# ╟─e7fd2c69-097d-400f-aa20-19265ebdd2eb
# ╠═987b5101-f2d6-4d03-98de-595bf6710d96
# ╟─a7cc3418-4607-4e03-ab87-bab26530cc53
# ╟─3651b884-9583-49fb-95e8-269f349cf1ce
# ╟─053b3718-7c29-4a81-9cb8-c639e515aa07
# ╟─f9e82af8-cbcb-4eab-8b66-227a5f34ccef
# ╟─78fc998d-a2c8-448b-90f5-d822e3513e8b
# ╟─e8278f0e-74a6-4cc3-8164-269391162662
# ╟─a8848e7f-dfd2-438a-956e-daa91ac7ecad
# ╟─8cfd5ff2-cbf6-469e-af5a-050ecf2f3198
# ╟─a69d2160-3a9e-4038-9d87-5119e9e05466
# ╟─95e1b0f1-fec8-4c35-a527-a5351ed38d92
# ╟─fdae435c-4d85-453e-a9ea-8383d31c5fe5
# ╟─90bd960f-794d-496f-bf28-3b52abba90cc
# ╟─62fb93e7-4acc-4466-b932-ce6cfdaf8d61
# ╟─30833702-0444-4d4e-86c6-be97ae1456e7
# ╟─a3ad9375-b2bf-4b3b-9f84-d18d41265640
# ╟─f50aff72-d978-4c2a-8683-8127a14a4ea9
# ╟─4dd54f9d-eb3f-49b4-9d56-41c5397ba001
# ╟─ee859bf9-4d20-46f5-9616-c932962cdbe2
# ╟─7d947cce-355f-4168-a76d-df5929d03be5
# ╟─9be1a640-c90f-4171-ab54-f6926dba25be
# ╟─3b64007d-762c-4bc5-8751-81ccb69ef376
# ╟─9e6ed715-6bf8-4e60-b77e-f4f8e2118f02
# ╟─b11a67a2-b544-4996-9816-82bfbcd70a18
# ╟─e69ad96e-7796-4db4-9ae1-049ad0971f9c
# ╟─505497ed-7060-48fe-ba2d-69a31413c267
# ╟─1b28028a-92df-4d7a-942f-11ab9f4d06a7
# ╟─ddb60135-0438-4381-8284-053c464ec506
# ╟─668930ff-00d1-46b6-97a3-b27f5f6628c1
# ╟─f35de386-da3a-4cbf-89ae-3049218531df
# ╟─9d2676de-fb2e-4a2c-8701-ee823bca5f71
# ╟─00594dfc-5a8b-4551-b300-bb52eee81e04
# ╟─a4e97910-11e5-4ad5-b9d2-7d5635e2990b
# ╟─22cb99e8-a5e4-4a16-8670-1ef1ef6f7b39
# ╟─b468696b-bb1b-44f6-8777-623b1b327b0b
# ╟─9ac51023-96ca-4304-a8b0-af36c3c8f60e
# ╟─a46e742e-9869-478c-b7a8-99267ceb9116
# ╟─b28b570c-44f3-49e9-9b94-eb6f2ed89bbf
# ╟─8cf9969b-c6ce-46af-9770-effd72bdf06c
# ╟─4ed1a6ca-0d21-4257-87db-f39fff0d208f
# ╟─54b12ace-4743-49a5-9e43-8580ee43ca6a
# ╟─12aa745e-9361-4af3-8c8b-7a2ffa83e874
# ╟─e3a20e0f-1524-479f-83e8-6fc2093e320b
# ╟─7ecbc555-7c10-4002-a662-b3de16611269
# ╟─9d6e52a1-d573-4b50-a584-5e517991e8ed
# ╟─e8756fef-22eb-48d9-b430-24b2af6dea4f
# ╟─d70b02b7-9697-428f-8eb7-e75f329da362
# ╟─9291e2ac-56a3-41a0-88b8-d3d8ecb7e819
# ╟─66b6ed61-d870-46d4-a5ad-cb6afec3a9dc
# ╟─ffb5767a-24fe-411f-abaa-baac67eaaa3d
# ╟─2cbceb55-e953-4035-a0d3-e4ec236038ab
# ╠═e4f0570b-0913-4f40-8a9b-8afb5ce7cbcd
# ╟─e57d0aa5-661f-43e8-885c-fd0aa1ab4cc6
# ╟─8ce9ee8c-cca0-4ff5-a5a0-14991987feb0
# ╟─e0ab7d80-17fc-478e-90a1-4f0922bfd728
# ╟─00e4aeb4-d28c-4de6-b298-a47a9d8ee3ab
# ╟─4e765c64-eac2-4654-badf-222601c888b7
# ╟─079acdc9-781f-43c1-bd93-12e40142a0af
# ╟─f91c04c2-4e85-4709-a76e-639929d54abd
# ╠═65193530-df0b-46f5-8653-5379a16ab779
# ╟─c638088a-c6fe-406e-8ef5-f3511319aef9
# ╟─fd29ef1c-d683-4a84-9595-04d0503d61ab
# ╟─59929002-a335-4368-ab7c-a9ad080b3e78
# ╟─70ad83fd-6df4-47b7-b5c0-0225c552d2e7
# ╟─261a5231-b450-4945-ba85-28c737e5f46f
# ╟─eec05ea5-cb8a-4979-a5c9-23eebd9afe3f
# ╟─3e3022c4-06ad-4878-aa53-79e274dd40ec
# ╟─cc596f86-0762-4878-9c32-aaf3981b6398
# ╟─972e9347-72c0-4a53-b20a-205086a29f8b
# ╟─123cec56-4f6d-445d-9f53-629161a1d487
# ╟─54e544f0-f1d0-467d-9a7a-b773053895df
# ╟─26847ae3-9b6d-46b4-84e0-225177041c15
# ╟─1c83d413-5ab1-40c2-a533-25a312889f5c
# ╟─06126200-59c5-4277-bb0c-f156ed90f51f
# ╟─1c4799ce-911b-4980-80fc-7c56c2b2a6ff
# ╠═c101ca1b-13af-4d0b-806b-14481f81b13e
# ╟─90297730-1473-4336-a884-d9441f3103a9
# ╟─e680937f-1fc5-4416-9991-8153bf604d64
# ╠═4f9efa7b-a9ea-4012-87e9-7d0cedb9be54
# ╟─e8f8c6c1-184c-4edf-a0b7-ce13a6eaf2a0
# ╟─cfbdaefd-bdcc-43e5-bdb4-87b37071c3c8
# ╟─e3deaaf2-9c9b-4693-b903-4c3ecfda1cd7
# ╟─e68ab39e-fdb5-4f39-8f18-c25f65e99c05
# ╟─ea5cb02b-49b4-4edf-89ab-d2fd6155e2a2
# ╟─c04f7b8a-f3ad-49d7-9610-48c5a7305649
# ╠═3d5e1274-386b-4511-8353-0188b2b65eb0
# ╟─611bfd31-5eab-41ee-b410-120739748a2a
# ╠═5d5adcae-2ad4-44b7-af71-a8ee276df8a3
# ╟─c3fb95a4-c9c0-4fc6-85d4-f6cf8a8258e4
# ╟─6a1d699c-fc16-472d-8395-a486890a089d
# ╟─61b0d517-94fc-4c6c-814b-8dfe51da6f1b
# ╟─5f768b64-c67e-4470-bf44-35bb539e1441
# ╟─a9aa5b73-bde3-4907-971d-696e26ada195
# ╟─47e27a13-4e60-40db-9ed7-5433e0230bd3
# ╟─c4a80b01-c1e3-4196-8408-3054d0aca71e
# ╟─25317d39-d107-41e4-a5e3-5ee11a0f8bd5
# ╟─d65a3a75-9e63-4087-912a-89e5e973a171
# ╟─1fd86b21-8457-4c1f-961c-81a105bda5a3
# ╟─bff8bd64-7a8c-4105-96b1-a3f8e8f0cd41
# ╟─c5f4f53e-6200-4a02-ada7-72e8ae60e38b
# ╟─cc599ac5-7000-45c9-bb5a-4529df717046
# ╟─4177166a-d638-4739-91e7-d0c42a80392d
# ╠═c2c69250-1f18-4c6d-b8b7-3de26ba0818a
# ╟─be8b89bb-e38f-423d-80b3-fa00fac7ad9a
# ╟─2482f975-0c63-45ab-8f9a-a2752809b192
# ╠═e45f6f80-7a93-4064-97b0-d6e7269a1b68
# ╟─755a1304-ec77-4ed8-806e-c287f13fdb89
# ╟─3115c2f6-8574-4d09-af83-db8e09731d07
# ╟─01b33b12-1dad-45e9-a191-8f9619b2fac5
# ╟─b39bcd45-dfc3-4ddb-95e8-7272ba591682
# ╟─eb60c5be-bad0-4c39-9842-3210dbaba8d3
# ╟─1e34cf47-8338-484e-b27a-fabfc7ef1b0b
# ╠═0c695eed-ebe3-4a26-9a32-7054bf1a51c3
# ╟─f0d5e130-d667-4f38-bfef-5aca795b880f
# ╠═92244fae-7546-4608-83ad-6dfe7425243c
# ╟─5fff03aa-5063-4bd0-bad1-6a3e568b25a1
# ╟─f50ef025-5faa-4a17-9307-4687bd3daac6
# ╟─4fde53d9-ae78-495e-9898-65af518d6b7f
# ╟─44e027ce-9f36-4ffb-bf81-b10905107770
# ╟─0cfec2e7-9166-45ad-b7cd-9585a1205e22
# ╟─2c65f0af-38dd-4d5b-868d-c62b6f2653f9
# ╟─9219bf2c-6e18-4a82-bf59-5050701bc404
# ╠═f24491c0-11a1-45b0-adac-cb9ce23f300f
# ╠═941815a4-e9e9-4d9d-b019-31da9ff3b176
# ╟─35608f80-9cf7-441b-a4f8-b5617899528d
# ╟─4b236b75-12dd-4638-b968-4923123af93d
# ╟─766a993a-b41e-4de3-b511-6fa6e4b092e7
# ╠═3141d33a-677b-47bc-ba45-fb906c984259
# ╟─3ade9674-04a8-4592-bb30-e807ffb0e17d
# ╠═aa7ba0e8-c1a3-4bdf-8bf5-aa389bc3795e
# ╟─92f9bc62-9d0f-40ee-832c-dff52a5d935b
# ╟─ac719200-b17f-4be9-be43-43b37b5c4018
# ╟─4171e8d8-5821-4a11-ba03-b8f88d785e68
# ╟─479e56b4-cee6-461e-8e87-21127f8272f5
# ╟─e30e4b18-9d86-46d1-89f3-c15c156aff6a
# ╟─61de967c-5ba7-44db-8aad-175026d4a62d
# ╟─15ba7f52-c201-4c7f-be84-da0bdf6a75f9
# ╟─34a1ae00-ae62-4782-bb03-9de7a8d26d1f
# ╠═4909326c-c521-40f6-9941-ade4f56adfa7
# ╠═9c2ec805-5efb-4b27-aa80-825220df4d66
# ╟─c7f882dd-faa0-45f5-bf2e-9100c72094ba
# ╟─619ac568-e272-4e55-8dfc-c965f4b878c1
# ╟─4a7f2e34-fe1e-4a34-a942-a4cf29938167
# ╟─bf7b6d08-4801-4757-87f0-c61e8966fec3
# ╟─ac71abdd-790e-403d-97f3-c88fe62b40b7
# ╟─4f4fcb2c-a636-40e3-b82e-a2b8888ca7ae
# ╟─b6163d7f-1447-461b-8912-f3c152fad7d5
# ╠═ce129375-fc38-46ae-a9b7-90605211620e
# ╠═c85a6afe-cb7f-4188-8950-d216956b7e7a
# ╟─e4e2c6aa-a849-4841-974e-a30e55843718
# ╟─c6947376-99f5-4564-822c-0578f136b54a
# ╟─f69d01bf-9832-4624-9c19-6c6eb024d24a
# ╠═997136cf-24c5-45b4-bec1-87d28504d7b0
# ╟─c1894e87-0df3-431f-8ce6-f4dab38349c7
# ╟─233e1c2d-2d45-4aa7-be03-6bc4dc76d716
# ╟─02e6a419-66f8-4aec-b456-6667b629308a
# ╟─979c7a6f-dd1d-4eb2-a862-7f370d654058
# ╟─0663cf1a-a42e-440f-931f-8321cc7d65f7
# ╟─e0c0c21a-8f1b-4b8e-8b4f-9fe8821c556b
# ╟─3119f567-fa6f-48a7-ac66-0cf630e3418f
# ╟─6e52e8cf-5b78-41bd-b8fd-390e37614b76
# ╟─e880c7b4-ab83-4bfc-ac7d-3d1616e3f925
# ╟─de6d1d86-94d3-4b8f-966c-4681dc4c5e2c
# ╟─d426f020-40fe-4f10-a366-5ec99c11900b
# ╟─00f5e3a3-dce4-4478-bf7c-2d1845205a71
# ╠═e7eed220-a52f-43d4-b431-57344875c0f7
# ╟─76715703-febf-454f-85b3-9bca7d8a44f4
# ╟─dda247d8-e3b5-4d0f-a072-4a6f0dc262b5
# ╟─32bf866c-03ec-48d4-bf41-31586e248a92
# ╟─8e77e17a-99a0-4c5d-8aa3-6e86f89f6292
# ╠═76bc523f-cbc4-49e0-9dbd-286b2a98a894
# ╟─b086377d-2be1-485f-a241-12df11cde33e
# ╟─4f1c3f5a-1490-4a1d-b26b-679773b90fac
# ╟─ab5d5e0d-afb3-457d-a169-e34608c1a93e
# ╟─17bcabbc-02c5-451e-9cbe-6e964cd306ef
# ╠═4e814e85-56f6-4407-9de8-73a4c0d30eb6
# ╟─07786494-0587-436d-a270-6f8e43d944e2
# ╟─0318bdb6-cec0-4f1f-bd0d-89c10553b211
# ╟─47d7efd0-48f6-479c-836a-f19c45d1318a
# ╟─ba9178c0-677b-470e-8164-7f54bd2438b5
# ╟─9ad7f65f-1e92-425c-9aa0-9804bbf1911d
# ╟─e60aaf53-b008-40aa-8403-1f619b81e5b7
# ╟─427b3198-6485-462b-8d13-23322c96a8ac
# ╟─a92e9f8f-06de-4176-aa1a-aed8a24bc71a
# ╟─6444c4a4-3015-46d1-8c22-8b28f54ad102
# ╟─0a5b7cc2-a84e-403e-bc7d-b40f94550383
# ╟─5189936d-6808-4d88-81da-454b992079eb
# ╟─4c56c9df-9186-4369-9021-c049a2eefdd0
# ╟─3d590d35-f866-404d-9c2a-e8dc17712113
# ╟─c6a776c3-c345-44d7-83f3-088f7eeb20b8
# ╟─74d5b989-04e0-45db-a5c5-316541e04fd7
# ╟─45df680f-d7a0-4112-9234-fdd69fa4f92c
# ╟─20ffd955-9b04-4562-b06b-8ee38a4abaa9
# ╟─cf37eeeb-7317-4050-8218-2a67a2432852
# ╠═842361ff-b144-463f-8ae1-9be8a2682723
# ╟─8566080c-a6bc-4823-b1f2-747687cdada1
# ╠═a554a842-3b8d-45b9-9d40-466bbb879986
# ╟─66820a81-3ad6-4905-9810-0ab5dd38ce69
# ╟─91cf0689-7e38-4c5d-b718-b444c73853f6
# ╠═e41760c1-8b60-4efa-9bce-dd915fd9b671
# ╟─b77f106b-5943-4053-9afc-a91a1554781b
# ╟─bf206834-929c-43fd-a35d-9cd1aa2976b2
# ╟─0c1b7f88-2907-474c-877b-4cfaf212dcf7
# ╟─a231aed6-741a-4fc7-b687-65042726dc3b
# ╟─fc5b69cd-eddd-4f1a-92c4-df4d26939821
# ╠═99eefc06-0c8a-4326-af71-45aa94814703
# ╟─008d89d6-7cdc-4d37-bf36-50efed0d03be
# ╟─c2169908-379d-4e5d-9f7f-d6fcecdf8f20
# ╠═767b00e9-3d6c-413c-8a82-1e9152aa5052
# ╠═65bdc140-2f92-11ef-1cbe-31065d820068
# ╟─4dd7bade-7523-4fa6-a862-25d2c61dbf9a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
