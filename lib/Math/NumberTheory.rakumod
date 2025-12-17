use v6.d;

unit module Math::NumberTheory;

#==========================================================
# Integer digits
#==========================================================
#| Gives a list of the base $b digits of the integer $n.
proto sub integer-digits(
        $n,                           #| An integer or list of integers to find the digits of.
        Int:D $base where 2..36 = 10, #| Base of the digits.
                   ) is export {*}

multi sub integer-digits(Int:D $n, Int:D $base = 10) {
    return $n.base($base).comb.map({ $_.Str.parse-base($base) }).List;
}

multi sub integer-digits(@n, Int:D $base = 10) {
    return @n.map({ integer-digits($_, $base) }).List;
}

#==========================================================
# GCD
#==========================================================
# It is a good idea to make GCD work for Gaussian integers and rationals.

proto sub gcd-gaussian(|) is export {*}

multi sub gcd-gaussian(*@n) {
    return reduce(&gcd-gaussian, @n);
}
multi sub gcd-gaussian(Int:D $a, Complex:D $b) {
    return gcd-gaussian($a + 0i, $b);
}
multi sub gcd-gaussian(Complex:D $a, Int:D $b) {
    return gcd-gaussian($a, $b + 0i);
}

multi sub gcd-gaussian(Complex:D $a is copy, Complex:D $b is copy) {
    sub norm(Complex:D $z) {
        return $z.re ** 2 + $z.im ** 2;
    }

    sub round-to-nearest-int($x) {
        return $x.round;
    }

    sub divide-gaussian(Complex:D $a, Complex:D $b) {
        my $denom = norm($b);
        my $real-part = ($a.re * $b.re + $a.im * $b.im) / $denom;
        my $imag-part = ($a.im * $b.re - $a.re * $b.im) / $denom;
        my $q-re = round-to-nearest-int($real-part);
        my $q-im = round-to-nearest-int($imag-part);
        return Complex.new($q-re, $q-im);
    }

    # Is this needed?
    if norm($b) > norm($a) {
        ($a, $b) = $b, $a;
    }

    while $b != 0 {
        my $q = divide-gaussian($a, $b);
        my $r = $a - $b * $q;
        $a = $b;
        $b = $r;
    }

    return $a;
}

#----------------------------------------------------------
proto sub gcd-rational(|) is export {*}

multi sub gcd-rational(*@r) {
    return reduce(&gcd-rational, @r);
}

multi sub gcd-rational(Rat $r1, Rat $r2) {
    my $num-gcd = $r1.numerator gcd $r2.numerator;
    my $den-lcm = $r1.denominator lcm $r2.denominator;
    return $num-gcd / $den-lcm;
}

#----------------------------------------------------------
# Redefine gcd

multi infix:<gcd>(Int:D $a, Complex:D $b --> Complex:D) is export {
    return gcd-gaussian($a, $b);
}

multi infix:<gcd>(Complex:D $a, Int:D $b --> Complex:D) is export {
    return gcd-gaussian($a, $b);
}

multi infix:<gcd>(Complex:D $a, Complex:D $b --> Complex:D) is export {
    return gcd-gaussian($a, $b);
}

multi infix:<gcd>(Rat:D $a, Rat:D $b --> Rat:D) is export {
    return gcd-rational($a, $b);
}

multi infix:<gcd>(Int:D $a, Rat:D $b --> Rat:D) is export {
    return gcd-rational($a.Rat, $b);
}

multi infix:<gcd>(Rat:D $a, Int:D $b --> Rat:D) is export {
    return gcd-rational($a, $b.Rat);
}

#==========================================================
# LCM
#==========================================================

proto sub lcm-gaussian(|) is export {*}

multi sub lcm-gaussian(*@n) {
    return reduce(&lcm-gaussian, @n);
}
multi sub lcm-gaussian(Int:D $a, Complex:D $b) {
    return lcm-gaussian($a + 0i, $b);
}
multi sub lcm-gaussian(Complex:D $a, Int:D $b) {
    return lcm-gaussian($a, $b + 0i);
}

multi sub lcm-gaussian(Complex:D $a is copy, Complex:D $b is copy) {
    return $a * $b / gcd-gaussian($a, $b);
}

#----------------------------------------------------------

proto sub lcm-rational(|) is export {*}

multi sub lcm-rational(*@r) {
    return reduce(&lcm-rational, @r);
}

multi sub lcm-rational(Rat $r1, Rat $r2)  {
    my $numerator-lcm = $r1.numerator lcm $r2.numerator;
    my $denominator-gcd = $r1.denominator gcd $r2.denominator;
    return $numerator-lcm / $denominator-gcd;
}

#----------------------------------------------------------
# Redefine lcm

multi infix:<lcm>(Int:D $a, Complex:D $b --> Complex:D) is export {
    return lcm-gaussian($a, $b);
}

multi infix:<lcm>(Complex:D $a, Int:D $b --> Complex:D) is export {
    return lcm-gaussian($a, $b);
}

multi infix:<lcm>(Complex:D $a, Complex:D $b --> Complex:D) is export {
    return lcm-gaussian($a, $b);
}

multi infix:<lcm>(Rat:D $a, Rat:D $b --> Rat:D) is export {
    return lcm-rational($a, $b);
}

multi infix:<lcm>(Int:D $a, Rat:D $b --> Rat:D) is export {
    return lcm-rational($a.Rat, $b);
}

multi infix:<lcm>(Rat:D $a, Int:D $b --> Rat:D) is export {
    return lcm-rational($a, $b.Rat);
}

#==========================================================
# PrimeQ
#==========================================================
# Extending is-prime to deal with Gaussian Integers.

#multi sub is-prime(Complex:D $p, Bool:D :gaussian(:$gaussian-integers) = True --> Bool) is export {...}
#| Gives True if the argument is Gaussian prime.
#| C<:$p> -- An integer or complex number or a list of numbers.
proto sub is-prime-gaussian($p --> Bool) is export {*}

multi sub is-prime-gaussian(Int:D $p --> Bool) {
    return is-prime($p) && $p mod 4 == 3;
}
multi sub is-prime-gaussian(Complex:D $p --> Bool) {
    return do if $p.re == 0 {
        is-prime($p.im.Int) && $p.im.Int mod 4 == 3
    } elsif $p.im == 0 {
        is-prime($p.re.Int) && $p.re.Int mod 4 == 3
    } else {
        my $a = $p.re.Int ** 2 + $p.im.Int ** 2;
        is-prime($a)
    }
}

multi sub is-prime-gaussian(@p --> List) {
    return @p.map({ is-prime-gaussian($_) }).List
}

#----------------------------------------------------------
# Redefine is-prime
multi sub is-prime(Complex:D $p) is export {
    return is-prime-gaussian($p);
}
multi sub is-prime($p, Bool:D :gaussian(:$gaussian-integers)) is export {
    return is-prime-gaussian($p);
}
multi sub is-prime(@p, *%args --> List) {
    return @p.map({ is-prime($_, |%args) }).List
}

#==========================================================
# Factorial
#==========================================================
# This has to be refactored -- in "Math::SpecialFunctions" has also &factorial.
# &factorial should work on reals:
# factorial($n) = gamma(1+$n)
#| Give the factorial of the argument.
#| C<:$n> -- Integer.
multi sub factorial(Int:D $n) is export {
    ([*] 1 .. $n) or 1
}

#==========================================================
# Fibonacci
#==========================================================
# This can be cached either by using "Memoize" or by
#  use experimental :cached
# and declaration as
#  sub fibonacci(UInt $n) is cached {...}
#| Give the n-th Fibonacci number.
#| C<$n> --  Non-negative integer.
sub fibonacci(UInt $n) is export {
    return 0 if $n == 0;
    return 1 if $n ≤ 2;
    my $k = $n div 2;
    if $n %% 2 {
        my $a = fibonacci($k);
        return $a * (2 * fibonacci($k + 1) - $a);
    }
    return fibonacci($k + 1) ** 2 + fibonacci($k) ** 2;
}

#==========================================================
# Digit count
#==========================================================
# http://reference.wolfram.com/language/ref/DigitCount.html
#| gives the number of digits in a base representation of the argument.
#| C<$n> -- Integer
#| C<:$base> - Base (integer)
#| C<:$digits> -
proto sub digit-count(Int:D $n, |) is export {*}
multi sub digit-count(Int:D $n, Int:D :b(:$base) = 10, :d(:$digits) = Whatever) {
    return do given $digits {
        when Whatever {
            my %h = $n.base($base).Str.comb.classify(*)».elems;
            %h
        }
        when Int:D {
            $n.base($base).Str.comb.grep($digits).elems
        }
        when $_ ~~ (List:D | Array:D | Seq:D) && $_.all ~~ Int:D {
            my $s = Set.new($digits».Str);
            my %h = $n.base($base).Str.comb.classify(*).map({ $_.key ∈ $s ?? ($_.key => $_.value.elems) !! Empty });
            %h
        }
        default {
            die 'The arugment $digit is expected to be Whatever, an integer, or a list of integers.'
        }
    }
}

multi sub digit-count(Int:D $n, Int:D $base = 10, $digits = Whatever) {
    return digit-count($n, :$base, :$digits)
}

#==========================================================
# Integer factors
#==========================================================

proto sub factor-integer($n, $k = Whatever, :$method = Whatever, Bool:D :gaussian(:$gaussian-integers) = False) is export {*}
multi sub factor-integer(
        $n is copy,
        $k is copy = Whatever,
        :$method is copy = Whatever,
        Bool:D :gaussian(:$gaussian-integers) = False) {
    die 'The first argument is expected to be an integer or a complex number with integer real and imaginary parts.'
    unless $n ~~ Int:D || ($n ~~ Complex:D && $n.re.floor == $n.re && $n.im.floor == $n.im);

    if $k.isa(Whatever) { $k = Inf }
    die 'The second argument is expected to be a positive integer or Whatever.'
    unless $k ~~ Int:D && $k > 0 || $k ~~ Inf;

    if $method.isa(Whatever) { $method = 'rho' }
    die 'The argument $methos is expected to be Whatever or one of "pollard-rho" or "trial-division".'
    unless $method ∈ <rho pollard pollard-rho trial trial-division>;

    if $gaussian-integers || $n ~~ Complex:D {
        if $n ~~ Int:D { $n = $n + 0i}
        return factor-gaussian-integer($n);
    } else {
        if $n < 0 {
            return [|(-1, 1), |factor-integer(abs($n), $k, :$method, :!gaussian-integers)].List;
        }

        return trial-factor-integer($n, $k) if $k < Inf;

        return do given $method {
            when $_.lc ∈ <rho pollard pollard-rho> {
                my @res = rho-prime-factors($n);
                if $k ≤ @res.elems {
                    @res.head($k)
                } else {
                    $k < Inf ?? @res.head($k) !! @res
                }
            }
            when $_.lc ∈ <trial trial-division> { trial-factor-integer($n, $k) }
        }
    }
}

sub factors-of(UInt:D $n is copy, UInt:D $p = 2) {
    my @factors;
    my $count = 0;
    while $n %% $p {
        $n div= $p;
        $count++;
    }
    @factors.push(($p, $count)».clone) if $count > 0;
    return $n, @factors;
}

#----------------------------------------------------------
sub trial-factor-integer(Int $n is copy, $k = Inf) is export {

    return ((1, 1),) if $n == 1;

    my @res = factors-of($n, 2);
    my @factors = |@res.tail;
    $n = @res.head;

    my $d = 3;
    while $d.Num ** 2 <= $n && @factors.elems < $k {
        if is-prime($d) {
            my $count = 0;
            while $n %% $d {
                $n div= $d;
                $count++;
            }
            @factors.push(($d, $count)».clone) if $count > 0;
        }
        $d++;
    }
    @factors.push(($n.clone, 1)) if $n > 1 && @factors.elems < $k;
    return @factors;
}

#----------------------------------------------------------
sub pollard-rho(Int $n, Int :$seed = 2, Int :$c = 1) {
    my $x = $seed;
    my $y = $seed;
    my $d = 1;
    while $d == 1 {
        $x = ($x * $x + $c) % $n;
        $y = ($y * $y + $c) % $n;
        $y = ($y * $y + $c) % $n;
        $d = ($x - $y).abs gcd $n;
    }
    return $d;
}

sub rho-factor-integer(UInt:D $n is copy, Int :$seed = 2, Int :$c = 1) {
    return ((1, 1),) if $n == 1;

    my @res = factors-of($n, 2);
    my @factors = |@res.tail;
    $n = @res.head;

    my &factorize = sub (Int $n is copy, UInt $c = 1) {
        return if $n == 1;
        if $n.is-prime {
            @factors.push(($n.clone, 1));
            return;
        }
        my $factor = pollard-rho($n, :$seed, :$c);
        my $exp = 0;
        while $n %% $factor {
            $n div= $factor;
            $exp++;
        }
        @factors.push(($factor, $exp)».clone);
        &factorize($n);
    }

    &factorize($n);
    return @factors.sort(*.head).List;
}

#----------------------------------------------------------
# From RosettaCode: https://rosettacode.org/wiki/Prime_decomposition#Pure_Raku
# With the comment:
# This is a pure Raku version that uses no outside libraries.
# It uses a variant of Pollard's rho factoring algorithm and is fairly performant when factoring numbers < 2⁸⁰;
# typically taking well under a second on an i7. It starts to slow down with larger numbers,
# but really bogs down factoring numbers that have more than 1 factor larger than about 2⁴⁰.

sub rho-prime-factors(Int:D $n) {
    return ((1, 1),) if $n == 1;
    my @res = prime-factors($n);
    @res = @res.classify(*).map({ ($_.key, $_.value.elems) }).sort(*.head);
    return @res;
}

sub prime-factors (Int $n where *> 0) {
    return $n if $n.is-prime;
    return () if $n == 1;
    my $factor = find-factor($n);
    sort flat ($factor, $n div $factor).map: &prime-factors;
}

sub find-factor (Int $n, $constant = 1) {
    return 2 unless $n +& 1;
    if (my $gcd = $n gcd 6541380665835015) > 1 { # magic number: [*] primes 3 .. 43
        return $gcd if $gcd != $n
    }
    my $x = 2;
    my $rho = 1;
    my $factor = 1;
    while $factor == 1 {
        $rho = $rho +< 1;
        my $fixed = $x;
        my int $i = 0;
        while $i < $rho {
            $x = ($x * $x + $constant) % $n;
            $factor = ($x - $fixed) gcd $n;
            last if 1 < $factor;
            $i = $i + 1;
        }
    }
    $factor = find-factor($n, $constant + 1) if $n == $factor;
    $factor;
}

#----------------------------------------------------------
# First implementation, kept for reference.
#`[
sub factor-gaussian-integer(Complex:D $n) {
    # See https://codegolf.stackexchange.com/a/185311
    my @res;
    sub f2($_){{$!=0+|sqrt .abs²-$^a²;{($!=$_/my \w=$^b+$a*i)==$!.floor&&.abs>w.abs>1>return f2 w&$!}for -$!..$!}for ^.abs;@res.push($_)};
    f2($n);
    @res = @res.classify(*).map({ ($_.key, $_.value.elems) }).sort(*.head.abs);
    return @res;
}
]

#----------------------------------------------------------
# Re-programming of the Python implementation given here:
#   https://github.com/johnhw/GaussianFactorisation/blob/master/gaussian_factorise.py
# Seems the to follow the description here:
#   https://stackoverflow.com/a/2271645

sub factor-gaussian-integer(Complex:D $a is copy, Bool:D :$include-unit = True) {
    sub norm(Complex:D $z) {
        return $z.re ** 2 + $z.im ** 2;
    }

    my $n = norm($a);
    my @factors = |factor-integer($n.Int).map({ $_.head xx $_.tail }).flat;

    my @z-factors;

    while @factors.elems > 0 {
        my $factor = @factors.shift;

        my $u;
        if $factor == 2 {
            # either 1+1j or 1-1j; 1+1j chosen here
            $u = 1 + 1i;
        } elsif $factor mod 4 == 3 {
            # x = 3 mod 4, remove two copies of this factor
            $u = $factor;
            # Remove repeated factor (note assumes factors are in order)!
            @factors.shift;
        } else {
            # x = 1 mod 4
            # Find k, such that k^2 = -1 mod factor = (factor-1) mod factor
            my $n = (2 .. ($factor - 1)).pick;
            while power-mod($n, ($factor - 1) div 2, $factor) != $factor - 1 {
                $n = (2 .. ($factor - 1)).pick;
            }
            my $k = power-mod($n, ($factor - 1) div 4, $factor);

            # Try dividing in k+1j
            my $trial-factor = gcd-gaussian($factor, $k + 1i);
            my $q = ($a / $trial-factor).round;

            # If exact, we have a factor
            if norm($a - $q * $trial-factor) < 1e-12 {
                $u = $trial-factor;
            } else {
                # Otherwise it is the conjugate
                $u = $trial-factor.conj;
            }
        }

        # Track the remaining number so we have the final unit factor
        $a = $a / $u;
        @z-factors.push( $u.round );
    }

    # We might have a factor of -1, 1j, or -1j -- append it if requested
    my @res = $include-unit ?? @z-factors.push( $a.round ) !! @z-factors;

    @res = @res.classify(*).map({ ($_.key, $_.value.elems) }).sort(*.head.abs);
    return @res;
}

#==========================================================
# Divisors
#==========================================================

#| Give a list of the integers that divide the argument.
sub divisors($n) is export {
    gather do {
        take 1;
        if $n > 2 {
            for 2 .. ($n div 2) -> $i {
                take $i if $n %% $i;
            }
        }
        take $n if $n != 1;
    }
}

#| Give the divisor function σ(exp, n).
proto divisor-sigma($n, |) is export {*}
multi sub divisor-sigma($exponent, $n) {
    divisor-sigma($n, :$exponent)
}

multi sub divisor-sigma($n, Int:D :e(:$exponent) = 1) {
    [+] divisors($n).map(-> $j { $j ** $exponent });
}

#==========================================================
# Euler Phi
#==========================================================
# http://reference.wolfram.com/language/ref/EulerPhi.html
# https://mathworld.wolfram.com/TotientFunction.html
# https://en.wikipedia.org/wiki/Euler%27s_totient_function
#| Gives the Euler totient function ϕ(n) that counts the positive integers up to a given integer n
#| that are relatively prime to n.
#| C<$n> -- Number.
proto sub euler-phi($n) is export {*}
multi sub euler-phi($n) {
    (^$n).grep({ $_ gcd $n == 1 }).elems
}

multi sub euler-phi(@ns) {
    @ns.map({ euler-phi($_) }).List
}

constant &totient is export(:ALL) = &euler-phi;

#==========================================================
# Integer exponent
#==========================================================
# http://reference.wolfram.com/language/ref/IntegerExponent.html
#| Gives the highest power of b that divides n.
#| C<$n> -- Integer.
#| C<$b> -- (Base) integer.
sub integer-exponent(Int:D $n is copy, UInt:D $b= 10) is export {
    my $k = 0;
    while ($n mod $b) == 0 {
        $n = $n div $b;
        $k++
    }
    return $k;
}

#==========================================================
# Power mod
#==========================================================
# http://reference.wolfram.com/language/ref/PowerMod.html
# $b ^ e mod m
proto sub power-mod(Int:D $b, $e, Int:D $m) is export {*}

multi sub power-mod(Int:D $b is copy, @e is copy, Int:D $m) {
    return @e.map({ power-mod($b, $_, $m) }).List;
}
multi sub power-mod(Int:D $b is copy, Int:D $e is copy, Int:D $m) {
    my $res = try expmod($b, $e, $m);
    if $! { return Nil }
    return $res;
}

sub complex-mod(Complex:D $a, Int:D $m) {
    return ($a.re mod $m) + ($a.im mod $m) * i;
}

multi sub power-mod(Complex:D $b is copy, Int:D $e is copy, Int:D $m) {
    die 'If the first argument is a complex number, then it is expected to be a Gaussian integer.'
    unless ($b.re == $b.re.floor) && ($b.im == $b.im.floor);

    if $m == 1 { return 0 }
    my $r = 1;
    $b = complex-mod($b, $m);
    while $e > 0 {
        if $e mod 2 == 1 {
            $r = complex-mod($r * $b, $m)
        }
        $b = complex-mod($b * $b, $m);
        $e = floor($e / 2);
    }
    return $r;
}

multi sub power-mod(@b, Int:D $e is copy, Int:D $m) {
    @b.map({ power-mod($_, $e, $m) }).List
}

#==========================================================
# Modular inverse
#==========================================================
# https://reference.wolfram.com/language/ref/ModularInverse.html
# Find k^-1 such that k * k^-1 mod n = 1.
proto sub modular-inverse($k, Int:D $n) is export {*}

multi sub modular-inverse(@k, Int:D $n) {
    @k.map({ modular-inverse($_, $n) }).List
}

multi sub modular-inverse(Int:D $k, Int:D $n) {
    return do if $n < 0 {
        die "Negative modulus is not implemented (yet.)"
    } elsif $k < 0 {
        power-mod($n - abs($k) mod $n, -1, $n)
    } else {
        power-mod($k, -1, $n)
    }
}

#==========================================================
# Chinese remainder
#==========================================================
# http://reference.wolfram.com/language/ref/ChineseRemainder.html
#| Give the smallest x with x>=0 that satisfies all the integer congruences x mod m_i == r_i mod m_i.
#| C<:@r> -- List of remainders.
#| C<:@m> -- List of modules
#| C<:$d> -- Number, lower limit of the result.
sub chinese-remainder(@r, @m, $d = 0) is export {
    my $n = @r.elems;
    my $M = [*] @m;
    my $x = 0;

    for ^$n -> $i {
        my $Mi = $M div @m[$i];
        #my $yi = (1..@m[$i]).first(* * $Mi % @m[$i] == 1);
        my $yi = modular-inverse($Mi, @m[$i]);
        without $yi { return Nil }
        $x += @r[$i] * $Mi * $yi;
    }

    $x = $x % $M;
    $x += $M while $x < $d;
    return $x;
}

#==========================================================
# Primitive root
#==========================================================
# http://reference.wolfram.com/language/ref/PrimitiveRoot.html
#| Give primitive root of the argument.
proto sub primitive-root($n, :$method = Whatever) is export {*}

multi sub primitive-root(Int:D $n, :$method = Whatever) {
    my $phi = euler-phi($n);
    my @factors = factor-integer($phi, :$method)».head;

    for 2 .. $n - 1 -> $g {
        next unless $g gcd $n == 1;
        if [&&] @factors.map({ ($g ** ($phi div $_)) mod $n != 1 }) {
            return $g;
        }
    }
    return Nil;
}

multi sub primitive-root(@n, :$method = Whatever) {
    die 'If the first argument is a list then all ellements are expected to be integers.'
    unless @n.all ~~ Int:D;

    @n.map({ primitive-root($_, :$method) }).List
}

# &primitive-root should be refactored to use a version of &primitive-root-list's code.

# http://reference.wolfram.com/language/ref/PrimitiveRootList.html
#| Give a list of primitive roots of the argument.
proto sub primitive-root-list($n, :$method = Whatever) is export {*}

multi sub primitive-root-list(Int:D $n, :$method = Whatever) {
    my $phi = euler-phi($n);
    my @factors = factor-integer($phi, :$method)».head;

    # Make a faster no-primitive root check based on:
    # A primitive root exists if and only if $n is 1, 2, 4, p^k or 2*p^k, where p is an odd prime and k > 0.

    my @res;
    for 2 .. $n - 1 -> $g {
        next unless $g gcd $n == 1;
        if [&&] @factors.map({ ($g ** ($phi div $_)) mod $n != 1 }) {
            @res.push($g);
        }
    }
    return @res;
}

multi sub primitive-root-list(@n, :$method = Whatever) {
    die 'If the first argument is a list then all ellements are expected to be integers.'
    unless @n.all ~~ Int:D;

    @n.map({ primitive-root-list($_, :$method) }).List
}

#==========================================================
# Multiplicative order
#==========================================================
# Is this supposed to work with Gaussian primes?
#| Give the multiplicative order of $^a ($k) modulo $^b ($n).
#|
proto sub multiplicative-order(Int:D $k is copy, Int:D $n where * > 0, | --> Int) is export {*}

multi sub multiplicative-order(Int:D $k is copy, Int:D $n where * > 0 --> Int) {
    # Handle edge cases
    return 1 if $n == 1;

    # Ensure k and n are coprime (gcd = 1)
    my $gcd = $k gcd $n;
    die 'The first and second arguments must be coprime.' if $gcd != 1;

    # Get the number of coprimes
    my $phi = euler-phi($n);

    # Get prime factorization of φ(n)
    my @phi-factors = factor-integer($phi);

    # Start with maximum possible order
    my $order = $phi;

    # For each prime factor of φ(n), find the smallest exponent needed
    for @phi-factors -> [$p, $e] {
        # Try removing each power of p from the order
        for 1..$e {
            my $test-order = $order div $p;
            if expmod($k, $test-order, $n) == 1 {
                $order = $test-order;
            } else {
                last;
            }
        }
    }

    return $order;
}

# General case: find smallest m where k^m ≡ r_i (mod n)
multi sub multiplicative-order(Int:D $k, Int:D $n where * > 0, *@r --> Int) {
    # Handle edge cases
    return 1 if $n == 1;

    # Ensure k and n are coprime (gcd = 1)
    my $gcd = $k gcd $n;
    die 'The first and second arguments must be coprime.' if $gcd != 1;

    # Normalize remainders and ensure they're valid
    my @remainders = @r.map({ $_ % $n }).unique;
    die 'The third argument, if given, is expected to be a list of integers.'
    unless @remainders.all ~~ Int:D;

    # Corner case to delegate
    if @remainders.elems == 1 && @remainders[0] == 1 {
        return multiplicative-order($k, $n);
    }

    my $m = 0;
    my $value = 1;
    my %seen = 1 => True;  # Track values to detect cycles

    loop {
        if @remainders.any == $value {
            return $m;
        }

        $m++;
        $value = ($value * $k) mod $n;

        # If we've seen this value before, we're in a cycle
        if %seen{$value}:exists {
            die "No solution exists for given remainders.";
        }
        %seen{$value} = True;

        # Safety limit (should be less than n)
        if $m >= $n {
            die "No solution exists for given remainders.";
        }
    }
}

#==========================================================
# Prime
#==========================================================
# http://reference.wolfram.com/language/ref/Prime.html
#| Give the n-th prime number.
#| C<$n> -- Prime number index or a list of indexes.
proto sub prime($n) is export {*}

multi sub prime(@n) {
    die 'All elements of a positional argument are expected to be integers.'
    unless @n.all ~~ Int:D;
    return @n».&prime;
}

multi sub prime(Int:D $n) {
    die "The argument is expected to be a positive integer."
    unless $n > 0;

    my $i = 0;
    my $candidate = 0;
    while $i < $n {
        $candidate++;
        if $candidate.is-prime { $i++ }
    }
    return $candidate;
}

sub primes-up-to(UInt:D $n) {
    my @res;
    my $i = 0;
    my $candidate = 0;
    while $i < $n {
        $candidate++;
        if $candidate.is-prime { $i++; @res.push($candidate) }
    }
    return @res;
}

#==========================================================
# Next prime
#==========================================================
# http://reference.wolfram.com/language/ref/NextPrime.html
#| Give the smallest prime number above the argument.
#| C<:$x> -- A number or a list of numbers.
#| C<:$k> -- An integer or a list of integers.
proto sub next-prime($x, |) is export {*}

multi sub next-prime(@xs, $k) {
    next-prime(@xs, :$k)
}

multi sub next-prime(@xs, :offset(:$k) = 1) {
    @xs.map({ next-prime($_, :$k) }).List
}

multi sub next-prime(Numeric:D $x, $k) {
    next-prime($x, :$k)
}

multi sub next-prime(Numeric:D $x, :offset(:@k)!) {
    return Nil unless @k.all ~~ Int:D;
    @k.map({ next-prime($x, k => $_)}).List
}

multi sub next-prime(Numeric:D $x, Int:D :offset(:$k) = 1) {
    given $x {
        when $_ !~~ Int:D { next-prime($x.floor, $k) }
        when $_ ~~ Int:D && $k == 1 {
            if $_ < -3 || $_ > 1 {
                integer-next-prime($_)
            } else {
                given $_ {
                    when -3 { -2 }
                    when -2 { 2 }
                    when -1 { 2 }
                    when 0 { 2 }
                    when 1 { 2 }
                }
            }
        }
        when is-prime($_) && $k == 0 { $_ }
        when $k > 0 {
            my $res = $x;
            for ^$k { $res = next-prime($res, k => 1) };
            $res
        }
        when $k < 0 {
            my $res = $x;
            for ^(-$k) { $res = previous-prime($res) };
            $res
        }
        default {
            Nil
        }
    }
}
sub integer-next-prime(Int:D $n) {
    my $candidate = $n.floor + 1 + $n mod 2;
    while !abs($candidate).is-prime {
        $candidate += 2;
    }
    return $candidate;
}

sub previous-prime(Numeric:D $x) {
    return -next-prime(-$x);
}

#==========================================================
# PrimePi
#==========================================================
#| Give the number of primes π(x) less than or equal to the argument.
proto sub prime-pi(Numeric:D $x, :$method = Whatever) is export {*}

multi sub prime-pi(Numeric:D $x, :$method is copy = Whatever) {
    if $method.isa(Whatever) { $method = 'legendre' }
    die 'The value of $method is expected to be a string or Whatever.'
    unless $method ~~ Str:D;

    return do given $method.lc {
        when 'sieve' { (1 .. $x.floor).grep(*.&is-prime).elems }
        when $_ ∈ <legendre legendre-formula> { legendre-formula($x) }
        default {
            die "Unknown method."
        }
    }
}

#----------------------------------------------------------
# Legendre's formula counts the number of positive integers
# less than or equal to a number x which are not divisible
# by any of the first a primes.
# https://mathworld.wolfram.com/LegendresFormula.html

# Speeding up computations
state %legendre-formula-cache;
# Using the recurrence relation
sub legendre-formula(Numeric:D $x) {
    sub phi($x, $a) {
        return $x if $a == 0;
        return %legendre-formula-cache{$x}{$a} //= phi($x, $a - 1) - phi(($x div prime($a)), $a - 1);
    }

    sub pi($n) {
        return 0 if $n < 2;
        my $a = pi(($n.sqrt).floor);
        return phi($n, $a) + $a - 1;
    }

    return pi($x.floor);
}

#==========================================================
# Prime omega
#==========================================================
#| Give the number of prime factors counting multiplicities Ω(n) in the argument.
multi sub prime-omega($x) is export {*}
multi sub prime-omega(@x) {
    return @x.map({ prime-omega($_) }).List;
}
multi sub prime-omega(Int:D $x) {
    return $x == 1 ?? 0 !! factor-integer($x)».tail.sum;
}

multi sub prime-omega(Complex:D $x) {
    return $x == 1 ?? 0 !! factor-gaussian-integer($x, :!include-unit)».tail.sum;
}

#==========================================================
# Prime nu
#==========================================================
#| Give the number of distinct primes in the argument.
multi sub prime-nu($x) is export {*}
multi sub prime-nu(@x) {
    return @x.map({ prime-nu($_) }).List;
}
multi sub prime-nu(Int:D $x) {
    return $x == 1 ?? 0 !! factor-integer($x).elems;
}

multi sub prime-nu(Complex:D $x) {
    return $x == 1 ?? 0 !! factor-gaussian-integer($x, :!include-unit).elems;
}

#==========================================================
# MoebiusMu function
#==========================================================
# http://reference.wolfram.com/language/ref/MoebiusMu.html

# Give the Möbius function µ(n).
proto sub moebius-mu($n, Bool:D :gaussian(:$gaussian-integers) = False) is export {*}

multi sub moebius-mu(@n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return @n.map({ moebius-mu($_, :$gaussian-integers) }).List;
}

multi sub moebius-mu(Int:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return moebius-mu($n + 0i) if $gaussian-integers;
    return 1 if $n == 1;
    my @exps = factor-integer($n)».tail;
    return @exps.max == 1 ?? (-1) ** @exps.elems !! 0;
}

multi sub moebius-mu(Complex:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    my @exps =  factor-gaussian-integer($n, :!include-unit)».tail;
    return @exps.max == 1 ?? (-1) ** @exps.elems !! 0;
}

#==========================================================
# SquareFreeQ function
#==========================================================
# http://reference.wolfram.com/language/ref/SquareFreeQ.html

# Give True if the first argument is a square-free number, and False otherwise.
proto sub is-square-free($n, Bool:D :gaussian(:$gaussian-integers) = False) is export {*}

multi sub is-square-free(@n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return @n.map({ is-square-free($_, :$gaussian-integers) }).List;
}

multi sub is-square-free(Int:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return moebius-mu($n, :$gaussian-integers) != 0 ?? True !! False;
}

multi sub is-square-free(Complex:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return moebius-mu($n) != 0 ?? True !! False;
}

#==========================================================
# LiouvilleLambda function
#==========================================================
# http://reference.wolfram.com/language/ref/LiouvilleLambda.html

# Gives the Liouville lambda function λ(n).
proto sub liouville-lambda($n, Bool:D :gaussian(:$gaussian-integers) = False) is export {*}

multi sub liouville-lambda(@n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return @n.map({ liouville-lambda($_, :$gaussian-integers) }).List;
}

multi sub liouville-lambda(Int:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return liouville-lambda($n + 0i) if $gaussian-integers;
    return (-1) ** prime-omega($n);
}

multi sub liouville-lambda(Complex:D $n, Bool:D :gaussian(:$gaussian-integers) = False) {
    return (-1) ** prime-omega($n);
}

#==========================================================
# CarmichaelLambda function
#==========================================================
# http://reference.wolfram.com/language/ref/CarmichaelLambda.html
# https://en.wikipedia.org/wiki/Carmichael_function#Recurrence_for_λ(n)

#| Give the Carmichael function λ(n).
#| C($n) -- An integer.
proto sub carmichael-lambda($n) is export {*}

multi sub carmichael-lambda(@n) {
    return @n.map({ carmichael-lambda($_) }).List;
}
multi sub carmichael-lambda(Int:D $n is copy) {
    return 0 if $n == 0;

    $n = abs($n);

    return 1 if $n == 1;

    my @factors = factor-integer(abs($n));

    my @lambdas = do for @factors -> $p {
        my $r = ($p.head - 1) * $p.head ** ($p.tail - 1);
        if $p.head == 2 && $p.tail ≥ 3 {
           $r /= 2
        }
        $r
    }
    return [lcm] @lambdas;
}

#==========================================================
# Kronecker delta
#==========================================================
#| Gives the Kronecker delta value, which is equal to 1 if all the arguments are equal, and 0 otherwise.
proto sub kronecker-delta(**@n) is export {*}

multi sub kronecker-delta() {
    return 1;
}
multi sub kronecker-delta(Numeric:D $n) {
    return $n == 0 ?? 1 !! 0;
}

multi sub kronecker-delta(+@n) {
    return kronecker-delta() if @n.elems == 0;
    return kronecker-delta(@n.head) if @n.elems == 1;
    return ([&&] @n.tail(*-1).map({ $_ == @n.head })) ?? 1 !! 0;
}

#==========================================================
# Simple extensions
#==========================================================
#| Gives true if the number is not prime.
proto sub is-composite($n) is export {*}
multi sub is-composite(Complex:D $n) { !is-prime-gaussian($n) }
multi sub is-composite(Int:D $n) { !is-prime($n) }

#| Gives true the number is a power of a prime.
sub is-prime-power(Int:D $n) is export {
    factor-integer($n).elems == 1;
}

#| Give Mangold lambda for an integer
sub mangold-lambda(Int:D $n) is export {
    my @res = factor-integer($n);
    return @res.elems == 1 ?? @res.head.head.log !! 0;
}

#==========================================================
# Coprime numbers
#==========================================================
#| Give true if the arguments are coprime.
proto sub are-coprime(|) is export {*}

multi sub are-coprime($n1, $n2) {
    ($n1 gcd $n2) == 1
}

multi sub are-coprime(*@ns) {
    die 'All elements are expected to be numbers.'
    unless @ns.all ~~ Numeric:D;

    for @ns.kv -> $i, $n1 {
        for @ns.kv -> $j, $n2 {
            if $i != $j {
                return False if $n1 gcd $n2 > 1
            }
        }
    }
    return True;
}

#==========================================================
# Twin, cousin, and sexy primes
#=========================================================

sub related-primes(UInt:D $n, UInt:D :$step = 2) is export {
    # This is slow, but it would be nice if it is made to run faster.
    #my @ps = (1..(2 * $n * log($n)).floor)».&prime;

    my @ps = primes-up-to(floor(2 * $n * log($n)));
    my @res = (@ps (&) (@ps >>+>> $step)).keys.sort;
    @res = @res[^(min($n, @res.elems))];
    return @res.map({ ($_ - $step, $_) }).List;
}

#| Get the first n-pairs of primes that differ by 2.
sub twin-primes(UInt:D $n) is export {
    related-primes($n, :2step);
}

#| Get the first n-pairs of primes that differ by 4.
sub cousin-primes(UInt:D $n) is export {
    related-primes($n, :4step);
}

#| Get the first n-pairs of primes that differ by 6.
sub sexy-primes(UInt:D $n) is export {
    related-primes($n, :6step);
}

#==========================================================
# Happy numbers
#==========================================================
#| Test whether an integer is a happy number.
proto sub is-happy-number(
        $n,                           #= An integer or a list of integers to check.
        Int:D $base where 2..36 = 10, #= Base to the digits.
        Int:D $p where $p > 0 = 2,    #= Power to raise the digits with.
                    ) is export {*}

multi sub is-happy-number(
        Int:D $n,
        Int:D $base where 2..36 = 10,
        Int:D $p where $p > 0 = 2) {
    my @out = [$n, ];
    my %seen;
    my $a;
    repeat {
        $a = integer-digits(@out.tail, $base).map({ $_ ** $p }).sum;
        @out.push($a);
        %seen{$a}++;
    } while %seen{$a} ≤ 1;

    return @out.tail == 1;
}

multi sub is-happy-number(@n, Int:D $base = 10, Int:D $p = 2) {
    return @n.map({ is-happy-number($_, $base, $p) }).List;
}

#==========================================================
# Integer partitions
#==========================================================
sub accel-asc(Int:D $n) {
    my @a = 0 xx ($n + 1);
    my $k = 1;
    my $y = $n - 1;

    gather {
        while $k != 0 {
            my $x = @a[$k - 1] + 1;
            $k--;

            while 2 * $x <= $y {
                @a[$k] = $x;
                $y -= $x;
                $k++;
            }

            my $l = $k + 1;
            while $x <= $y {
                @a[$k] = $x;
                @a[$l] = $y;
                take @a[0..($k + 1)];
                $x++;
                $y--;
            }

            @a[$k] = $x + $y;
            $y = $x + $y - 1;
            take @a[0..$k];
        }
    }
}

#| Give a list of all possible ways to partition the integer argument into smaller integers.
proto sub integer-paritions(Int:D $n, |) is export {*}

multi sub integer-paritions(Int:D $n, UInt:D $k) {
    integer-paritions($n, k-max => $k);
}

multi sub integer-paritions(Int:D $n, (Numeric:D $k-min, Numeric:D $k-max)) {
    integer-paritions($n, :$k-min, :$k-max);
}

multi sub integer-paritions(Int:D $n, UInt:D :$k-min = 1, Numeric:D :$k-max = Inf) {
    my @res = accel-asc($n);
    @res .= grep({ $k-min ≤ $_.elems ≤ $k-max });
    return @res».reverse».List.reverse.List;
}

#==========================================================
# Random prime
#==========================================================
# http://reference.wolfram.com/language/ref/RandomPrime.html
#| Give pseudo random prime numbers in specified ranges.
#| C<:$range> -- Integer or Range.
#| C<:$n> -- Number of primes to return.
proto sub random-prime(|) is export {*}
multi sub random-prime(Int:D $max, $n = Whatever) {
    die 'If the first argument is a number then it is expected to be an integer greater than 1.'
    unless $max > 1;
    return random-prime(2 .. $max, $n);
}

multi sub random-prime(Range:D $range, $n is copy = Whatever) {
    my ($min, $max) = $range.head, $range.tail;

    die 'Correct range argument is expected.'
    unless $min.defined && $max.defined;

    die 'The start of the range argument is expected to start with an integer greater than 1.'
    unless $min > 1;

    die 'The end of the range argument is expected to be greater than the range start.'
    unless $min < $max;

    die 'The second argument is expected to be a positive integer or Whatever.'
    unless $n ~~ Int:D && $n > 0 || $n.isa(Whatever);

    my $i = 0;
    my @res;
    my $n2 = $n.isa(Whatever) ?? 1 !! $n;
    for $min .. $max -> $c {
        if $c.is-prime {
            @res.push($c);
            $i++
        }
        last if $i ≥ $n2
    }

    return do if $n.isa(Whatever) {
        @res.head
    } elsif @res.elems == $n2 {
        @res.List
    } else {
        @res.roll($n2).List
    }
}

#==========================================================
# Real digits
#==========================================================
# http://reference.wolfram.com/language/ref/RealDigits.html
#| Real digits
#| C<$x> -- Number to convert.
#| C<:$base> -- Conversion base.
#| C<:$n> -- Digit exponent to start with.
#| C<:$tol> -- Tolerance to stop the conversion with.
#| C<:$length> -- Max number of digits.
proto sub real-digits(Numeric:D $x, *@args, *%args) is export {*}

multi sub real-digits($x, Numeric:D $b = 10, $n = Whatever, *%args) {
    return real-digits($x ~~ Int:D ?? $x.FatRat !! $x, :$b, :$n, |%args);
}

multi sub real-digits(Numeric:D $x is copy,
                      Numeric:D :base(:$b) = 10,
                      :start(:first-digit-exponent(:$n)) is copy = Whatever,
                      Numeric:D :tolerance(:$tol) = 10e-14,
                      :l(:len(:$length)) is copy = Inf) {
    if $x == 0 { return 0, 0 }
    $x = abs($x);
    my @digits;
    my $current = $x;

    if $length.isa(Whatever) { $length = Inf }
    die 'The argument $length is expected to be a positive integer or Whatever.'
    unless ($length === Inf || $length ~~ Int:D) && $length > 0;

    if $n.isa(Whatever) { $n = $x.log($b).floor }
    die 'The third argument is expected to be a number or Whatever.' unless $n ~~ Numeric:D;

    if $x ~~ FatRat:D {
        my $exp = $n.FatRat;
        my $bf = $b.FatRat;

        my $bfe = [*] ($bf xx $exp);
        while $current / $x > $tol && @digits.elems < $length {
            my $r = $current / $bfe;
            my $digit = do if $r.round(10 ** -100) == 1 { 1 } else { ($current / $bfe).floor };
            @digits .= push($digit);
            $current -= $digit.FatRat * $bfe;
            $exp--;
            $bfe = $bfe / $bf;
        }
    } else {
        my $exp = $n;
        my $bf = $b;

        my $bfe = [*] ($bf xx $exp);
        while $current / $x > $tol && @digits.elems < $length {
            my $digit = floor($current / $bfe);
            @digits .= push($digit);
            $current -= $digit * $bfe;
            $exp--;
            $bfe = $bfe / $bf;
        }
    }

    return @digits, $n + 1;
}

#==========================================================
# Phi number system
#==========================================================
# https://resources.wolframcloud.com/FunctionRepository/resources/PhiNumberSystem/
# https://mathworld.wolfram.com/PhiNumberSystem.html

# Using N[Sqrt[5], 100] in Wolfram Language
constant $sqrt5 = 2.236067977499789696409173668731276235440618359611525724270897245410520925637804899414414408378782275.FatRat;

# Fibonacci 401 and 400
constant $fibonacci401 = 284812298108489611757988937681460995615380088782304890986477195645969271404032323901.FatRat;
constant $fibonacci400 = 176023680645013966468226945392411250770384383304492191886725992896575345044216019675.FatRat;

#| Golden ratio (phi)
our constant \phi is export = $fibonacci401 / $fibonacci400;
our constant \ϕ is export = $fibonacci401 / $fibonacci400;
#our constant \ϕ is export = (1.FatRat + $sqrt5) / 2.FatRat;
sub golden-ratio(Bool:D :pre(:$pre-computed) = True) is export {
    return $pre-computed ?? \phi !! (1 + sqrt(5)) / 2;
}

#| Phi number system.
#| C<$n> -- An integer number to convert.
#| C<:$tol> -- Tolerance to stop the conversion with.
#| C<:$length> -- Max number of digits.
sub phi-number-system(Int:D $n, Numeric:D :tolerance(:$tol) = 10e-16, :l(:len(:$length)) is copy = Whatever) is export {
    if $length.isa(Whatever) { $length = 2 * ($sqrt5 * abs($n)).log(ϕ).floor + 1; }
    my ($digits, $exp) = real-digits($n.FatRat, ϕ, :$tol, :$length);
    return $exp <<->> ($digits.grep(*== 1, :k) >>+>> 1);
}
