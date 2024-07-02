struct FloatK
    m::Int
    n::Int
    s::Int
    digits::Int
    rounding::String
end
Base.convert(::Type{Float64}, x::FloatK) = x.s * x.m * (10.0)^(x.n - x.digits)
Base.show(io::IO, ::MIME"text/plain", n::FloatK) = print(io, n.s == 1 ? "" : "-", " 0.", n.m, " × 10^", n.n)
Base.show(io::IO, ::MIME"application/x-tex", n::FloatK) = print(io, n.s == 1 ? "" : "-", "0.", n.m, "× 10^", n.n)
FloatK(num::T, digits::S, truncation::String) where {T<:Real,S<:Int} = begin
    s = sign(num)
    num_str = "$(Float64(abs(rationalize((num)))))"
    num_str_with_dot = if 'e' in num_str
        e_pos = findfirst(x -> x == 'e', num_str)
        power = parse(Int, num_str[e_pos+1:end])
        is_power_positive = power >= 0
        value = parse(Int, replace(num_str[1:e_pos-1], '.' => ""))
        num_zeros = is_power_positive ? abs(power) : abs(power) - 1
        zeros_str = repeat('0', num_zeros)
        is_power_positive ? "$value" * zeros_str * ".0" : "0." * zeros_str * "$value"
    else
        num_str
    end
    num_parts = split(num_str_with_dot, '.')
    integral_part, fraction_part = num_parts[1], num_parts[2]
    positive_exp, integral_part_str = if length(integral_part) == 1 && integral_part[1] == '0'
        0, ""
    else
        length(integral_part), integral_part
    end
    negative_exp, fraction_part_str = if length(fraction_part) == 1 && fraction_part[1] == '0'
        0, ""
    else
        -count(x -> x == '0', fraction_part), fraction_part
    end
    floatk_power = positive_exp + negative_exp
    floatk_power_str = "$(integral_part_str)$(fraction_part_str)" * repeat('0', digits)
    mantissa = if truncation == "chop"
        parse(Int, floatk_power_str[1:digits])
    else
        after_last_digit = length(floatk_power_str) > digits ? parse(Int, floatk_power_str[digits+1]) : 0
        correction = after_last_digit >= 5 ? 1 : 0
        last_digit = parse(Int, floatk_power_str[digits]) + correction
        parse(Int, floatk_power_str[1:digits-1] * "$(last_digit)")
    end

    # floatk_power, mantissa, positive_exp, negative_exp, num_str_with_dot
    FloatK(mantissa, floatk_power, s, digits, truncation)
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
Base.:-(a::FloatK, b::FloatK) = FloatK(Float64(rationalize(convert(Float64, a)) - rationalize(convert(Float64, b))), a.digits, a.rounding)
Base.:-(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, a) - b, a.digits, a.rounding)
Base.:-(b::T where {T<:Real}, a::FloatK) = FloatK(b - convert(Float64, a), a.digits, a.rounding)

Base.:*(a::FloatK, b::FloatK) = FloatK(Float64(rationalize(convert(Float64, a)) * rationalize(convert(Float64, b))), a.digits, a.rounding)
Base.:*(a::T where {T<:Real}, b::FloatK) = FloatK(a * convert(Float64, b), b.digits, b.rounding)
Base.:*(a::FloatK, b::T where {T<:Real}) = FloatK(b, a)

Base.:^(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, a)^b, a.digits, a.rounding)
Base.:sqrt(a::FloatK) = FloatK(sqrt(convert(Float64, a)), a.digits, a.rounding)

Base.:÷(a::FloatK, b::FloatK) = FloatK(Float64(rationalize(convert(Float64, a)) / rationalize(convert(Float64, b))), a.digits, a.rounding)
Base.:÷(a::T where {T<:Real}, b::FloatK) = FloatK(a / convert(Float64, b), b.digits, b.rounding)
Base.:÷(a::FloatK, b::T where {T<:Real}) = FloatK(convert(Float64, b) / a, a.digits, a.rounding)

Base.:/(a::FloatK, b::FloatK) = FloatK(Float64(rationalize(convert(Float64, a)) / rationalize(convert(Float64, b))), a.digits, a.rounding)
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

fl, ⊕, ⊖, ⊗, ⨸ = createFiniteDigitSystem(; digits=5, truncation="chop")
x = 5.0 / 7.0
y = 1.0 / 3.0
x ⊕ y
x ⊖ y
x ⊗ y
x ⨸ y

b = fl(1.1011)
f(x) = x^3 - 6.1x^2 + 3.2 * x + 1.5
fp(x) = fl(x)^3 - 6.1 * x^2 + 3.2 * fl(x) + 1.5
x = 4.71
xp = fl(x)
f(xp)
fp(x)
-6 * (xp * xp)
fl(0.0001)

x = 5 / 7
y = 1 / 3
u = 0.714251
v = 98765.9
w = 0.111111e-4
p1 = fl(x) - fl(u)
p2 = p1 ÷ fl(w)

function tofloatk(num::T, digits::S, truncation::String) where {T<:Real,S<:Int}
    s = sign(num)
    num_str = "$(Float64(abs(num)))"
    num_str_with_dot = if 'e' in num_str
        e_pos = findfirst(x -> x == 'e', num_str)
        power = parse(Int, num_str[e_pos+1:end])
        println(power)
        is_power_positive = power >= 0
        value = parse(Int, replace(num_str[1:e_pos-1], '.' => ""))
        num_zeros = is_power_positive ? abs(power) : abs(power) - 1
        zeros_str = repeat('0', num_zeros)
        is_power_positive ? "$value" * zeros_str * ".0" : "0." * zeros_str * "$value"
    else
        num_str
    end
    num_parts = split(num_str_with_dot, '.')
    integral_part, fraction_part = num_parts[1], num_parts[2]
    positive_exp, integral_part_str = if length(integral_part) == 1 && integral_part[1] == '0'
        0, ""
    else
        length(integral_part), integral_part
    end
    negative_exp, fraction_part_str = if length(fraction_part) == 1 && fraction_part[1] == '0'
        0, ""
    else
        -count(x -> x == '0', fraction_part), fraction_part
    end
    floatk_power = positive_exp + negative_exp
    floatk_power_str = "$(integral_part_str)$(fraction_part_str)" * repeat('0', digits)
    mantissa = if truncation == "chop"
        parse(Int, floatk_power_str[1:digits])
    else
        after_last_digit = length(floatk_power_str) > digits ? parse(Int, floatk_power_str[digits+1]) : 0
        correction = after_last_digit >= 5 ? 1 : 0
        last_digit = parse(Int, floatk_power_str[digits]) + correction
        parse(Int, floatk_power_str[1:digits-1] * "$(last_digit)")
    end

    floatk_power, mantissa, positive_exp, negative_exp, num_str_with_dot, num_str
    # FloatK(mantissa, floatk_power, s, digits, truncation)
end
fl(x)
fl(u)
xmu1 = rationalize(0.71428)
xmu2 = rationalize(0.71425)
xmu = xmu1 - xmu2
Float64(xmu)
Float64(rationalize(0.1 + 0.2))