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
	using LinearAlgebra, Random,  Printf
	using Symbolics
	using QRCoders
	using PrettyTables
	# using NonlinearSolve
	# using ForwardDiff
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
	
	exportqrcode("https://mshahrani.website/teaching/","website.png"; width=1)
	
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
	D(f,x) = begin
		expand_derivatives(Df(f(x)))
	end
	df(f,n,x0) = begin 
		val = if n== 0
		f(x)
		else
		reduce((acc,val)->begin 
			D(t->substitute(acc,Dict(x=>t)),x)
				end,2:n;init=D(f,x))
		end
		substitute(val,Dict(x=>x0))
	end
end

# ╔═╡ b54035ab-2813-4cdc-887a-e11625ede4aa
slider1h = @bind slider1 Slider(0:20, show_value=true);""

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
	x =split("27.56640625",'.')
	# sum(float(Int(x[2][i])//(2^(i))) for i in 1:length(x[2]))
	bs=bitstring(parse(Float64,"27.56640625"))
	reinterpret(Float64,parse(Int,"0100000000111011100100001111111111111111111111111111111111111111",base=2))
	
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
	aerror1=0.1e2
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

	function re(p::T where {T<:Number},ps::FloatK)
		abs(convert(Float64,(p-ps))/p)
	end
end

# ╔═╡ e223e89c-af0c-44ed-bb83-fd51124c7899
let
	f(x)=cos(x)
	P2(x)=1-0.5*x^2
	x=3
	f(x),P2(x)
end

# ╔═╡ ee13d4c2-edce-4d7c-8ace-103aab1f7ba0
let
	f(x) = cos(x)
	x0=0
	# slider1
	
	P(n)=x0->sum(((x-x0)^i)*df(f,i,x0)/factorial(i) for i in 0:n)
	Pn = P(slider1)(x0) 
	Pn	
	
end

# ╔═╡ afbd4725-9c76-4453-a62a-ebac185655a2
let
	@syms ζ::Real
	x0=0
	n=slider1
	f(x) = cos(x)
	((x-x0)^(n+1))*df(f,n+1,ζ)/factorial(n+1)
end

# ╔═╡ 422a3ccc-1e62-407b-bd3e-0703f1a5a33d
let
	f(x)=cos(x)
	x0=0
	P(n)=x0->sum(((x-x0)^i)*df(f,i,x0)/factorial(i) for i in 0:n)
	Pn = P(slider1)(x0) 
	p1= plot(f;framestyle=:origin,label=L"f(x)",line=(1,:red))
	labels = [L"P_{%$i}(x)" for i in 0:2:slider1]
	plot(p1,[t->substitute(P(i)(x0),Dict(x=>t)) for i in 0:2:slider1];
	label=reshape(labels,1,length(labels)),
		xlimits=(-12,12),ylimits=(-1,1.5)
	)
end

# ╔═╡ 9b0cd30c-3162-4423-b32f-88c9a01c1cbd
let
	x1=0.01
	n=slider1
	f(x)=cos(x)
	x0=0
	P(n)=x0->sum(((x-x0)^i)*df(f,i,x0)/factorial(i) for i in 0:n)
	vals = substitute(P(n)(x0),Dict(x=>x1)), f(x1)
	Rn=((x-x0)^(n+1))/factorial(n+1) |> y -> substitute(y,Dict(x=>x1))
	pntex = texeq("""
	P_{$n}($x1) = $(vals[1])
		""")
	ftex=texeq("f($x1) = $(vals[2])") 
	error=texeq("|\\textrm{Error}| \\leq $(Rn)")
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
		if length(str)>64
			error("The binary string is too long. It should be less than 64.")
		end
		if length(filter(x->x in ['0','1'],str))!=length(str)
			error("The string is not valid.")
		end
		reduce((c,v)->v[2]=='0' ? c : c+2^(Float64(v[1]-1)),zip(length(str):-1:1,str);init=0)
	end
	function float64_to_dec(str::String)
		if length(str)>64
			error("The binary string is too long. It should be less than 64.")
		end
		if length(filter(x->x in ['0','1'],str))!=length(str)
			error("The string is not valid.")
		end
		lstr =length(str) 
		str =  lstr < 63 ? repeat('0',64-length(str))*str : str
		# m,c,s= str[13:64], str[2:12], str[1]
		# mv = reduce((c,v)->v[2]=='0' ? c : c+1/2^(BigFloat(v[1])),zip(1:length(m),m);init=0)
		# cv = bint64_to_dec(c)
		# sv = s=='0' ? 1 : -1
		# # reduce((c,v)->v[2]=='0' ? c : c+2^(v[1]-1),zip(length(str):-1:1,str);init=0)
		
		# sv*2^(float(cv-1023))*(1+mv)
		reinterpret(Float64,parse(Int,str,base=2))
	end
end

# ╔═╡ 726a42a1-9746-4083-88f3-0661cd7396b5
bint64_to_dec("11111111111")-1023

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
	next_small="0100000000111011100100001111111111111111111111111111111111111111"
	next_largest="0100000000111011100100010000000000000000000000000000000000000001"
	pn = float64_to_dec(next_small)
	n=float64_to_dec(binstr)
	nn=float64_to_dec.(next_largest)
	# 0.5(nn-n)
end

# ╔═╡ bb4ed82e-3b7e-4898-987f-9a9bac49c205
let 	
	smallest_number=String(repeat('0',64))
	# smll_normalized =String(repeat('0',11)*repeat('0',52))
	# bitstring(1-1023)
	largest_number='0'*repeat('1',62)*'0'
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
strs = join(map(x->"$(x*5*10^(-4.0))",[0.1  0.5  100  1000  5000  9990 10000]),",")
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
	re(p::T,ps::S) where {T<:Number,S<:Number} = abs((p-ps)/p)
	fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=5, truncation="chop")
	x = 5/7
	y=1/3
	p,ps=x+y,	x ⊕ y
	re(p,convert(Float64,ps))
	# x ⊖ y
	# x ⊗ y
	# x ⨸ y
end

# ╔═╡ b50497e4-e0f9-4da6-bba3-48628cee0173
let
	re(p::T,ps::S) where {T<:Number,S<:Number} = abs((p-ps)/p)
	fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=5, truncation="chop")
	x=5/7
	y=1/3
	u=0.714251
	v=98765.9
	w=0.111111e-4
	p1 = 0.0003
	p2 = p1 ⊗ v
end

# ╔═╡ 379608bf-906d-45db-9cab-018a50fcbd07
let
	f(x)=x^2 +62.10x+1
	x1(a,b,c) =((-b+sqrt(b^2-4*a*c))/2a)
	x2(a,b,c) =((-b-sqrt(b^2-4*a*c))/2a)
	a,b,c =1, 62.10, 1
	x1(a,b,c)
	x2(a,b,c)
	aa=cc=1
	ba=62.10
	ab2= 0.3856e4
	# ab2= ba^2
	a4ac = 4
	asqr=0.6206e2
	# asqr=sqrt(ab2-a4ac)
	anum1=-0.0400
	ax1=anum1/2.0
	anum2=-0.1242e3
	ax2=anum2/2.0
	
	
end

# ╔═╡ 1922d625-ba5e-4ea2-b050-a3c64077c742
let
	fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(digits=5, truncation="round")
	f(x)=x^2 +62.10x+1
	x1(a,b,c) =((-b+sqrt(b^2-4*a*c))/2a)
	x1f(a,b,c) =((-fl(b)+fl(sqrt(fl(fl(fl(b)^2)-fl(4*fl(a)*fl(c))))))/fl(2a))
	x2(a,b,c) =((-b-sqrt(b^2-4*a*c))/2a)
	x2f(a,b,c) =((-fl(b)-fl(sqrt(fl(fl(fl(b)^2)-fl(4*fl(a)*fl(c))))))/fl(2a))
	a,b,c=1.0,62.10,1.0
	x1(a,b,c)
	x2f(a,b,c)
	# x2(a,b,c)
	# x2(fl(a),fl(b),fl(c))
end

# ╔═╡ 5ed7d51c-b5fa-4953-b2b4-c2040edced33
md"## Nested Arithmetic"

# ╔═╡ 987b5101-f2d6-4d03-98de-595bf6710d96
let
	fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=3, truncation="chop")
	f(x)=x^3-6.1*x^2+3.2*x+1.5
	fn(x)=x*(x*(x-6.1)+3.2)+1.5
	x=4.71
	y=f(x)
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
	P(n)=sum(iseven(i) ? (-1)^i*h^(i)/(factorial(i)) : 0 for i in 0:n) + (iseven(n) ? sin(ζ)*h^(n+1)/factorial(n+1) : cos(ζ)*h^(n+1)/factorial(n+1))
	p3 = P(3)
	
end

# ╔═╡ 4dd54f9d-eb3f-49b4-9d56-41c5397ba001
md"# 2.1 The Bisection Method"

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

# ╔═╡ 3b64007d-762c-4bc5-8751-81ccb69ef376
let
	f(x)=x^2*(x+4)-10
	
	# p,T = bisect(1,2,1e-4,14)
	# p,T
	# pretty_table(HTML,T[1:13,:], header=["n","an", "bn", "pn", "f(pn)"])
end

# ╔═╡ 9e6ed715-6bf8-4e60-b77e-f4f8e2118f02
let
	function bisect(a,b,TOL,N0)
		i=1
		FA=f(a)
		T=Matrix{Number}(undef,N0,5)
		while i<=N0
			p = a+(b-a)/2
			FP=f(p)
			T[i,:]=vcat(i,a,b,p,FP)
			if FP==0 || ((b-a)/2)<TOL
				TT = T[1:i,:]
				return p,TT
			end
			i=i+1
			if FA*FP>0
				a = p
				FA=FP
			else
				b=p
			end
		end
		@error("Maximum number of iterations reached")
	end
end

# ╔═╡ 505497ed-7060-48fe-ba2d-69a31413c267
md"# 2.2 Fixed-Point Iteration"

# ╔═╡ 1b28028a-92df-4d7a-942f-11ab9f4d06a7
begin 
	function fixed_point(g,x0,ϵ;maxiters=50)
		n = maxiters
		# r(n,x0)=reduce((c,_)->abs(c-g(c)) <ϵ ? c : g(c),1:n,init=x0)
		xs = Vector{Float64}(undef,n)
		xs[1]=x0
		# fixit(x)=g(x)
		i = 2
		while i <= n
			xs[i]= try 
				g(xs[i-1])
			catch
				NaN
			end
			if abs(1-xs[i-1]/xs[i]) < ϵ || isnan(xs[i])
				break
			end
			i += 1
		end
		# gx = [r(i,x0) for i in 1:n]
		last_i = min(i,n)
		xs = filter(x->!isnan(x),xs[1:last_i])
		ys = xs[2:end]
		xss = xs[1:end-1] 
		xss,ys
	end
	function animate_fixedpoint(g,x0,ϵ;maxiters=50)
		xs,ys = fixed_point(g,x0,ϵ;maxiters=50)
		
		plt= plot([g,x->x],label=[L"g(x)" L"y=x"], framestyle=:origin,legend=:topright)
		anim = @animate for j ∈ 1:(length(xs)-1)
			scatter(plt,xs[1:j],ys[1:j],label=L"x_{%$j}=%$(xs[j])")
		end
		# # annotation=[(2,5,L"x_{%$i}=%$(xs[i])",10)]
		# # gif(anim, "anim_fps15.gif", fps = 2)
		anim
	end
	
end

# ╔═╡ f35de386-da3a-4cbf-89ae-3049218531df
# let
# 	anim = animate_fixedpoint(x->x^2-2,1.1,0.001)
# 	gif(anim, "anim1_fps15.gif", fps = 2)
# end

# ╔═╡ 9ac51023-96ca-4304-a8b0-af36c3c8f60e
md"## Fixed-Point Iteration"

# ╔═╡ 4ed1a6ca-0d21-4257-87db-f39fff0d208f
let
	p0=1.5
	n = 30
	ϵ = 1e-9
	g1(x)=x-x^3-4x^2+10
	g2(x)=sqrt((10/x)-4x)
	g3(x)=(1/2.0)*sqrt(10-x^3)
	g4(x)=sqrt(10/(4+x))
	g5(x)=x-(x^3+4x^2-10)/(3x^2+8x)
	xs1,ys1 = fixed_point(g1,p0,ϵ;maxiters=n)
	xs2,ys2 = fixed_point(g2,p0,ϵ;maxiters=n)
	xs3,ys3 = fixed_point(g3,p0,ϵ;maxiters=n)
	xs4,ys4 = fixed_point(g4,p0,ϵ;maxiters=n)
	xs5,ys5 = fixed_point(g5,p0,ϵ;maxiters=n)
	
	T = hcat(0:n,
		vcat(xs1,repeat([" "],n+1-length(xs1))),
		vcat(xs2,repeat([" "],n+1-length(xs2))),
		vcat(xs3,repeat([" "],n+1-length(xs3))),
		vcat(xs4,repeat([" "],n+1-length(xs4))),
		vcat(xs5,repeat([" "],n+1-length(xs5))))
	# T = hcat(reshape(repeat([" "],31),31,1))
	# length(xs1)
	pretty_table(HTML,T,header=vcat(:n, [Symbol("g$i") for i in 1:5]...))
	# anim2 = animate_fixedpoint(g4,p0,ϵ;maxiters=n)
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

# ╔═╡ 4dd7bade-7523-4fa6-a862-25d2c61dbf9a
begin
	function post_img(img::String,w=500)
		res=Resource(img,:width=>w)
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
		beginBlock("Definition",t)
	end
	function bbl(t)
		beginBlock(t,"")
	end
	function bbl(t,s)
		beginBlock(t,s)
	end
	ebl()=endBlock()
	function bth(s)
		beginTheorem(s)
	end
	eth()=endTheorem()
	ex(n::Int;s::String="")=ex("Example $n",s)
	ex(t,s)=example(t,s)
	function beginBlock(title,subtitle)
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
		beginBlock("Theorem",subtitle)
	end
	function endBlock()
		"""</p></div></div>"""
	end
	function endTheorem()
		 endBlock()
	end
	function example(lable,desc)
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
	text_book= post_img("https://m.media-amazon.com/images/I/51ziKPbuEmL.jpg",200);
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
Let ``g \in C[a, b]`` be such that ``g(x) \in[a, b]``, for all ``x`` in ``[a, b]``. Suppose, in addition, that ``g^{\prime}`` exists on ``(a, b)`` and that a constant ``0<k<1`` exists with
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QRCoders = "f42e9828-16f3-11ed-2883-9126170b272d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Colors = "~0.12.11"
CommonMark = "~0.8.12"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Latexify = "~0.16.3"
PlotThemes = "~3.1.0"
Plots = "~1.40.4"
PlutoExtras = "~0.7.12"
PlutoUI = "~0.7.59"
PrettyTables = "~2.3.2"
QRCoders = "~1.4.5"
Symbolics = "~5.28.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "2d3e45879b83a196527e2219b673a6adc424622e"

[[deps.ADTypes]]
git-tree-sha1 = "fa0822e5baee6e23081c2685ae27265dabee23d8"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.4.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

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

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

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

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

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

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

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

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

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
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

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

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

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
git-tree-sha1 = "20339c0dd70abdb73494955df4fcd9e9ccaff861"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.6.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "ForwardDiff", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils"]
git-tree-sha1 = "4104548fff14d7370b278ee767651d6ec61eb195"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.28.0"

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
git-tree-sha1 = "a947ea21087caba0a798c5e494d0bb78e3a1a3a0"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.9"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

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
# ╟─7d947cce-355f-4168-a76d-df5929d03be5
# ╟─9be1a640-c90f-4171-ab54-f6926dba25be
# ╠═3b64007d-762c-4bc5-8751-81ccb69ef376
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
# ╟─9ac51023-96ca-4304-a8b0-af36c3c8f60e
# ╟─a46e742e-9869-478c-b7a8-99267ceb9116
# ╟─b28b570c-44f3-49e9-9b94-eb6f2ed89bbf
# ╟─4ed1a6ca-0d21-4257-87db-f39fff0d208f
# ╟─54b12ace-4743-49a5-9e43-8580ee43ca6a
# ╟─12aa745e-9361-4af3-8c8b-7a2ffa83e874
# ╟─e3a20e0f-1524-479f-83e8-6fc2093e320b
# ╟─7ecbc555-7c10-4002-a662-b3de16611269
# ╟─9d6e52a1-d573-4b50-a584-5e517991e8ed
# ╠═65bdc140-2f92-11ef-1cbe-31065d820068
# ╟─4dd7bade-7523-4fa6-a862-25d2c61dbf9a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
