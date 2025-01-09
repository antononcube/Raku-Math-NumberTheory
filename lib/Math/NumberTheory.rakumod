use v6.d;

unit module Math::NumberTheory;

#==========================================================
# Popular support functions
#==========================================================

sub factorial($n) is export { ([*] 1..$n) or 1 }

#==========================================================
# Integer factors
#==========================================================

sub integer-factors(Int:D $n is copy, UInt:D $b= 10) is export {

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

sub divisor-sigma($n, $exponent=1) is export {
    [+] divisors($n).map( -> $j { $j ** $exponent });
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
sub power-mod(Int:D $b is copy, Int:D $e is copy, Int:D $m) is export {
    if $m == 1 {
        return 0
    }
    my $r = 1;
    $b = $b mod $m;
    while $e > 0 {
        if $e mod 2 == 1 {
            $r = ($r * $b) % $m
        }
        $b = ($b * $b) % $m;
        $e = floor($e / 2);
    }
    return $r;
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
#| C<$n> -- an integer number to convert.
#| C<:$tol> -- tolerance to stop the conversion with.
#| C<:$length> -- max number of digits.
sub phi-number-system(Int:D $n, Numeric:D :tolerance(:$tol) = 10e-16, :l(:len(:$length)) is copy = Whatever) is export {
    if $length.isa(Whatever) { $length = 2 * ($sqrt5 * abs($n)).log(ϕ).floor + 1; }
    my ($digits, $exp) = real-digits($n.FatRat, ϕ, :$tol, :$length);
    return $exp <<->> ($digits.grep(*== 1, :k) >>+>> 1);
}