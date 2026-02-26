use v6.d;

unit module Math::NumberTheory::ContinuedFraction;


#==========================================================
# From Generalized Continued Fraction
#==========================================================

#| Reconstruct a number from the list of its generalized continued fraction terms.
our proto sub from-generalized-continued-fraction(@a, @b) is export {*}
multi sub from-generalized-continued-fraction(@a, @b, :$n is copy = Whatever) {
    die 'The named argument $n is expected to be a non-negative integer or Whatever.'
    unless $n.isa(Whatever) || $n ~~ Int:D && $n ≥ 0;

    $n = do if @a.is-lazy && @b.is-lazy {
        die 'When both positional arguments are lazy then the named argument $n is expected to be a non-negative integer.'
        unless $n ~~ Int:D && $n ≥ 0;
        $n
    } elsif !@a.is-lazy && @b.is-lazy {
        $n.isa(Whatever) ?? @a.elems - 1 !! min(@a.elems - 1, $n)
    } elsif @a.is-lazy && !@b.is-lazy {
        $n.isa(Whatever) ?? @b.elems !! min(@b.elems, $n)
    } else {
        $n.isa(Whatever) ?? min(@a.elems - 1, @b.elems) !! min(@a.elems - 1, @b.elems, $n)
    }

    return @a.head if $n == 0;

    die 'The two positional arguments are expected to be non-empty positionals.' if $n < 0;

    note "{$n - @a.elems + 1} tailing elements of the first positional argument are ignored."
    if !@a.is-lazy && @a.elems - 1 < $n;

    note "{$n - @b.elems} tailing elements of the second positional argument are ignored."
    if !@b.is-lazy && @b.elems < $n;

    die 'The arguments are expected to be sequences of integers'
    unless @a.all ~~ Int:D && @b.all ~~ Int:D;

    my $x = @a[$n - 1];
    $x = @a[$_ - 1] + @b[$_] / $x for reverse 1 ..^ $n;
    return $x;
}

# Taken from https://rosettacode.org/wiki/Continued_fraction#Raku
#multi sub from-continued-fraction(@a, @b) {
#    die 'The arguments are expected to be sequences of integers'
#    unless @a.all ~~ Int:D && @b.all ~~ Int:D;
#
#    return map { .(Inf) }, [\o] map { @a[$_] + @b[$_] / * }, ^Inf;
#}

#==========================================================
# From Continued Fraction
#==========================================================
# http://reference.wolfram.com/language/ref/ContinuedFraction.html
# https://mathworld.wolfram.com/ContinuedFraction.html

#| Reconstruct a number from the list of its continued fraction terms.
our proto sub from-continued-fraction(@a, |) is export {*}

multi sub from-continued-fraction(@a) {
    die 'The first argument is expected to be sequences of integers.'
    unless @a.all ~~ Int:D;
    return from-generalized-continued-fraction(@a, (1 xx (@a.elems - 1)));
}

#==========================================================
# To Continued Fraction
#==========================================================
#| Generate a list of the first $n terms in the continued fraction representation of $x.
our proto sub continued-fraction($x, | --> List) is export {*}

multi sub continued-fraction($x, $n, :tol(:$tolerance) is copy = Whatever --> List) {
    return continued-fraction($x, :$n, :$tolerance)
}

multi sub continued-fraction(@x, :$number-of-terms = Whatever, :$tolerance = Whatever) {
    return @x.map({ continued-fraction($_, :$number-of-terms, :$tolerance) }).List;
}

multi sub continued-fraction(
        Numeric:D $x,
        :n(:$number-of-terms) is copy = Whatever,
        :tol(:$tolerance) is copy = Whatever
        --> List) {

    if $x < 0 {
        return -1 <<*>> continued-fraction(-$x, :$number-of-terms, :$tolerance);
    }

    die "The number of terms argument is expected to be a positive integer or Whatever."
    unless $number-of-terms ~~ Int:D && $number-of-terms > 0 || $number-of-terms.isa(Whatever);

    die "The tolerance argument is expected to be a non-negative number or Whatever."
    unless $tolerance ~~ Numeric:D && $tolerance ≥ 0 || $tolerance.isa(Whatever);

    if $number-of-terms.isa(Whatever) && $tolerance.isa(Whatever) {
        $number-of-terms = Inf;
        $tolerance = 10e-12;
    } elsif $tolerance.isa(Whatever) {
        $tolerance = 0;
    } elsif $number-of-terms.isa(Whatever) {
        $number-of-terms = Inf
    }

    # This uses the wrong "Cauchy sequence" tolerance check.
    # Has to be replaced with convergents denominators based check.
    my $r = $x;
    my $k = 0;
    my @res;
    loop {
        my $a = $r.floor;
        @res.push($a);
        my $f = $r - $a;
        last if $f == 0 || $k++ ≥ $number-of-terms - 1 || $tolerance && abs($x - from-continued-fraction(@res)) ≤ $tolerance;
        $r = 1 / $f;
    }
    return @res.List;
}

#==========================================================
# Convergents
#==========================================================
# http://reference.wolfram.com/language/ref/Convergents.html

#| Give a list of the convergents corresponding to the given continued fraction terms list or a number.
our proto sub convergents($x, | --> List) is export {*}

multi sub convergents($x, $n, :tol(:$tolerance) is copy = Whatever --> List) {
    return convergents($x, :$n, :$tolerance)
}

multi sub convergents(@x, :n(:$number-of-terms) = Whatever, :tol(:$tolerance) = Whatever --> List) {

    return @x.List if @x.elems < 2;

    # Rat or FatRat?
    my @res = @x[0].Rat, (1 + @x[0] * @x[1]) / @x[1];
    for @x[2..^*] -> $d {
        my $p = $d * @res[*-1].numerator + @res[*-2].numerator;
        my $q = $d * @res[*-1].denominator + @res[*-2].denominator;
        @res.push($p / $q);
    }
    return @res.List;
}

multi sub convergents(
        Numeric:D $x,
        :n(:$number-of-terms) is copy = Whatever,
        :tol(:$tolerance) is copy = Whatever
        --> List) {
    my @terms = continued-fraction($x, :$number-of-terms, :$tolerance);
    return convergents(@terms);
}
