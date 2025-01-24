use v6.d;

unit module Math::NumberTheory;

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

sub factor-integer(Int:D $n is copy, $k is copy = Whatever, :$method is copy = Whatever) is export {
    if $n < 0 {
        return [|(-1, 1), |factor-integer(abs($n), $k, :$method)].List;
    }

    if $k.isa(Whatever) { $k = Inf }
    die 'The second argument is expected to be a positive integer or Whatever.'
    unless $k ~~ Int:D && $k > 0 || $k ~~ Inf;

    if $method.isa(Whatever) { $method = 'rho' }
    die 'The argument $methos is expected to be Whatever or one of "pollard-rho" or "trial-division".'
    unless $method ∈ <rho pollard pollard-rho trial trial-division>;

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

#==========================================================
# Divisors
#==========================================================

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

multi sub power-mod(Complex:D $b is copy, Int:D $e is copy, Int:D $m) is export {
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

#==========================================================
# Modular inverse
#==========================================================
# https://reference.wolfram.com/language/ref/ModularInverse.html
# Find k^-1 such that k * k^-1 mod n = 1.
sub modular-inverse(Int:D $k, Int:D $n) is export {
    return do if $n < 0 {
        die "Negative modulus is not implemented (yet.)"
    } elsif $k < 0 {
        power-mod($n - abs($k) mod $n, -1, $n)
    } else {
        power-mod($k, -1, $n)
    }
}

#==========================================================
# Chinese reminder
#==========================================================
# http://reference.wolfram.com/language/ref/ChineseRemainder.html
#| Give the smallest x with x>=0 that satisfies all the integer congruences x mod m_i == r_i mod m_i.
#| C<:@r> -- List of reminders.
#| C<:@m> -- List of modules
#| C<:$d> -- Number, lower limit of the result.
sub chinese-reminder(@r, @m, $d = 0) is export {
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
sub primitive-root(Int:D $n, :$method = Whatever) is export {
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

# http://reference.wolfram.com/language/ref/PrimitiveRootList.html
#| Give a list of primitive roots of the argument.
sub primitive-root-list(Int:D $n, :$method = Whatever) is export {
    my $phi = euler-phi($n);
    my @factors = factor-integer($phi, :$method)».head;

    my @res;
    for 2 .. $n - 1 -> $g {
        next unless $g gcd $n == 1;
        if [&&] @factors.map({ ($g ** ($phi div $_)) mod $n != 1 }) {
            @res.push($g);
        }
    }
    return @res;
}

#==========================================================
# Prime
#==========================================================
# http://reference.wolfram.com/language/ref/Prime.html
#| Give the n-th prime number.
#| C<$n> -- Prime number index.
sub prime(Int:D $n) is export {
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

#==========================================================
# Next prime
#==========================================================
# http://reference.wolfram.com/language/ref/NextPrime.html
#| Give the smallest prime number above the argument.
#| C<:$x> -- Number.
proto sub next-prime(Numeric:D $x, |) is export {*}
multi sub next-prime(Numeric:D $x, Int:D $k) {
    next-prime($x, :$k)
}

multi sub next-prime(@xs, Int:D :$k = 1) {
    @xs.map({ next-prime($_, :$k) }).List
}

multi sub next-prime(Numeric:D $x, Int:D :$k = 1) {
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
    return factor-integer($x)».tail.sum;
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
    return factor-integer($x).elems;
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
constant $sqrt5 = 2.236067977499789696409173668731276235440618359611525724270897245410520925637804899414414408378782275
        .FatRat;

# Fibonacci 401 and 400
constant $fibonacci401 = 284812298108489611757988937681460995615380088782304890986477195645969271404032323901.FatRat;
constant $fibonacci400 = 176023680645013966468226945392411250770384383304492191886725992896575345044216019675.FatRat;

#| Golden ratio (phi)
our constant \phi is export = $fibonacci401 / $fibonacci400;
our constant \ϕ is export = $fibonacci401 / $fibonacci400;
#our constant \ϕ is export = (1.FatRat + $sqrt5) / 2.FatRat;

#| Phi number system.
#| C<$n> -- An integer number to convert.
#| C<:$tol> -- Tolerance to stop the conversion with.
#| C<:$length> -- Max number of digits.
sub phi-number-system(Int:D $n, Numeric:D :tolerance(:$tol) = 10e-16, :l(:len(:$length)) is copy = Whatever) is export {
    if $length.isa(Whatever) { $length = 2 * ($sqrt5 * abs($n)).log(ϕ).floor + 1; }
    my ($digits, $exp) = real-digits($n.FatRat, ϕ, :$tol, :$length);
    return $exp <<->> ($digits.grep(*== 1, :k) >>+>> 1);
}