use v6.d;

unit module Math::NumberTheory;

#==========================================================
# Popular support functions
#==========================================================

sub factorial($n) is export {
    ([*] 1 .. $n) or 1
}

#==========================================================
# Integer factors
#==========================================================

sub factor-integer(Int $n is copy, $k is copy = Whatever, :$method is copy = Whatever) is export {
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
            my @res = rho-factor-integer($n);
            if $k ≤ @res.elems {
                @res.head($k)
            } else {
                @res = @res.map(-> @f { rho-factor-integer(@f.head, seed => 3, c => 17).map({ ($_.head, $_.tail * @f.tail) }) }).map(*.Slip);
                @res = @res.classify(*.head).map({ ($_.key, $_.value.map(*.tail).sum) }).sort(*.head);
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

sub trial-factor-integer(Int $n is copy, $k is copy = Whatever) is export {

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
    @factors.push(($n.clone, 1)) if $n > 1;
    return @factors;
}

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

sub divisor-sigma($n, $exponent= 1) is export {
    [+] divisors($n).map(-> $j { $j ** $exponent });
}

#==========================================================
# Euler Phi
#==========================================================
# http://reference.wolfram.com/language/ref/EulerPhi.html
sub euler-phi ($n) is export() {
    +(^$n).grep({ $_ gcd $n == 1 })
}

#==========================================================
# Integer exponent
#==========================================================
# http://reference.wolfram.com/language/ref/IntegerExponent.html
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
    } elsif $k < 0  {
        power-mod($n - abs($k) mod $n, -1, $n)
    } else {
        power-mod($k, -1, $n)
    }
}

#==========================================================
# Chinese reminder
#==========================================================

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
# Prime
#==========================================================

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

sub next-prime(Numeric:D $x) is export {
    my $candidate = $x.floor + 1;
    while !abs($candidate).is-prime {
        $candidate++;
    }
    return $candidate;
}

#==========================================================
# Next prime
#==========================================================

proto sub random-prime(|) is export {*}
multi sub random-prime(Int:D $max, $n = Whatever) {
    die 'If the first argument is a number then it is expected to be an integer greater than 1.'
    unless $max > 1;
    return random-prime(2..$max, $n);
}

multi sub random-prime(Range:D $range, $n is copy = Whatever) {

    my ($min, $max) = $range.head, $range.tail;

    die 'The start of the range argument is expected to start with an integer greater than 1.'
    unless $max > 1;

    die 'The end of the range argument is expected to be greater than the range start.'
    unless $min < $max ;

    die 'The second argument is expected to be a positive integer or Whatever.'
    unless $n ~~ Int:D && $n > 0 || $n === Inf || $n.isa(Whatever);

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
#| C<$x> -- number to convert.
#| C<:$base> -- conversion base.
#| C<:$n> -- digit exponent to start with.
#| C<:$tol> -- tolerance to stop the conversion with.
#| C<:$length> -- max number of digits.
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
            #note (:$exp, power => $bf ** $exp, :$bfe);
            #note '$current / $bfe => ', $current / $bfe;
            my $r = $current / $bfe;
            my $digit = do if $r.round(10 ** -100) == 1 { 1 } else { ($current / $bfe).floor };
            @digits .= push($digit);
            #note (:$digit, '$digit.FatRat * $bfe => ', $digit.FatRat * $bfe);
            $current -= $digit.FatRat * $bfe;
            #note (:$current);
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

#| Phi number system.
#| C<$n> -- an integer number to convert.
#| C<:$tol> -- tolerance to stop the conversion with.
#| C<:$length> -- max number of digits.
sub phi-number-system(Int:D $n, Numeric:D :tolerance(:$tol) = 10e-16, :l(:len(:$length)) is copy = Whatever) is export {
    if $length.isa(Whatever) { $length = 2 * ($sqrt5 * abs($n)).log(ϕ).floor + 1; }
    my ($digits, $exp) = real-digits($n.FatRat, ϕ, :$tol, :$length);
    return $exp <<->> ($digits.grep(*== 1, :k) >>+>> 1);
}