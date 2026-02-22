unit module Math::NumberTheory::Expansions;

use Math::NumberTheory;

=begin pod
=head1 NAME

Math::NumberTheory::Expansion

=head1 SYNOPSIS

    use Math::NumberTheory::Expansion :ALL;

    say number-expansion(pi, t => "Engel", :10n);                        # → [1, 1, 1, 8, ...]
    say number-expansion(3/19, type => "Sylvester");                     # → [7, 67]
    say number-expansion(384, "Cantor");                                 # → [0, 0, 0, 1, 3]
    say number-expansion(pi, t => "CantorProduct", :5number-of-terms);   # → [1, 2, 22, 600, ...]
    say number-expansion(5/13, "Lueroth");                               # → [3, [3,4,2]]  (periodic)

    say from-number-expansion([2, 7, 3, 2, 2, 2, 2, 2], 'lueroth');      # approx. reconstruction
    say from-number-expansion([1,1,1,8,8], type => "Engel");             # approx. reconstruction
    say from-number-expansion((1, 0, 1, 0, 1, 0, 1), 'Zeckendorf');      # exact reconstruction
=end pod

# -------------------------
# Public API
# -------------------------

proto sub number-expansion($x, | --> List) is export {*}
multi sub number-expansion($x, Str:D $type, $n = Whatever --> List) {
    return do given $n {
        when $_.isa(Whatever) { number-expansion-all($x, $type) }
        when $_ ~~ Int:D && $_ ≥ 0 { number-expansion-n($x, $type, $n) }
        default {
            die 'The third positional argument is expected to be a non-negative integer or Whatever.'
        }
    }
}
multi sub number-expansion($x, Str:D :t(:$type), :n(:$number-of-terms) --> List) { number-expansion-n($x, $type, $number-of-terms) }
multi sub number-expansion($x, Str:D :t(:$type) --> List) { number-expansion-all($x, $type) }

proto sub from-number-expansion(@terms, |) is export {*}

multi sub from-number-expansion(@terms, Str:D $type) {
    return from-number-expansion(@terms, :$type);
}

multi sub from-number-expansion(@terms, Str:D :t(:$type) is copy) {
    $type = normalize-type($type);
    given $type {
        when 'Engel'         { return from-engel(@terms) }
        when 'Sylvester'     { return from-sylvester(@terms) }
        when 'Cantor'        { return from-cantor(@terms) }
        when 'CantorProduct' { return from-cantor-product(@terms) }
        when 'Lueroth'       { return from-lueroth(@terms) }
        when 'Zeckendorf'    { return from-zeckendorf(@terms) }
    default {
        die "from-number-expansion: unsupported expansion type '$type'.";
    }
}
}

# -------------------------
# Normalization / guards
# -------------------------

sub normalize-type(Str:D $t --> Str) {
    my $u = $t.trim.lc;
    return 'Engel'         if $u eq 'engel';
    return 'Pierce'        if $u eq 'pierce';
    return 'Sylvester'     if $u eq 'sylvester';
    return 'Cantor'        if $u eq 'cantor';
    return 'CantorProduct' if $u ∈ <cantorproduct cantor-product>;
    return 'Lueroth'       if $u ∈ <lueroth lüroth luroth lüroth>;
    return 'Zeckendorf'    if $u eq 'zeckendorf';
    return 'Oppenheim'     if $u eq 'oppenheim';
    $t.trim
}

sub is-exact-rational($x --> Bool) {
    # Rat / FatRat / Rational-ish; also Int counts as rational.
    ($x ~~ Int) || ($x ~~ Rat) || ($x ~~ FatRat)
}

sub require-finite-for-irrational($x, Str $type) {
    die "number-expansion($type): expanding an irrational number requires specifying a finite length n"
    unless is-exact-rational($x);
}

sub as-fat-rat($x --> FatRat) {
    # Best-effort exactification.
    return $x.FatRat if $x ~~ Rat | FatRat | Int;
    die "Expected an exact rational (Rat/FatRat/Int), got {$x.^name}";
}

# -------------------------
# Dispatcher
# -------------------------

sub number-expansion-n($x, Str:D $t, Int:D $n --> List) {
    die "n must be >= 0" if $n < 0;
    my $type = normalize-type($t);

    given $type {
        when 'Engel'         { return take-n(engel-terms($x),         $n) }
        when 'Pierce'        { return take-n(pierce-terms($x),        $n) }
        when 'Sylvester'     { return take-n(sylvester-terms($x),     $n) }
        when 'Cantor'        { return cantor-expansion($x) } # ignores n; finite digits
        when 'CantorProduct' { return take-n(cantor-product-terms($x),$n) }
        when 'Lueroth'       { return take-n(lueroth-terms($x),       $n) }
        when 'Zeckendorf'    { return zeckendorf($x) } # ignores n; finite bits
        when 'Oppenheim'     { die "Use number-expansion(x, n, r, s, p, q, \"Oppenheim\") for Oppenheim expansion." }
        default              { die "Unknown expansion type '$t' (normalized to '$type')" }
    }
}

sub number-expansion-all($x, Str:D $t --> Mu) {
    my $type = normalize-type($t);

    given $type {
        when 'Cantor' {
            die "Cantor expansion expects a non-negative Int" unless $x ~~ Int && $x >= 0;
            return cantor-expansion($x);
        }
        when 'Zeckendorf' {
            die "Zeckendorf representation expects a non-negative Int" unless $x ~~ Int && $x >= 0;
            return zeckendorf($x);
        }
        when 'Lueroth' {
            require-finite-for-irrational($x, $type);
            return lueroth-all-rational(as-fat-rat($x));
        }
        default {
            # Spec: for exact numbers, number-expansion(x,t) can be used if x is rational;
            # irrationals must specify n.
            require-finite-for-irrational($x, $type);
            my $fr = as-fat-rat($x);

            my $seq = do given $type {
                when 'Engel'         { engel-terms($fr) }
                when 'Pierce'        { pierce-terms($fr) }
                when 'Sylvester'     { sylvester-terms($fr) }
                when 'CantorProduct' { cantor-product-terms($fr) }
                default {
                    die "number-expansion(x, \"$t\"): type '$type' requires an explicit length n or is unsupported in all-terms form";
                }
            };

            # Consume until termination (exact rationals should terminate for these).
            return $seq.List;
        }
    }
}

sub take-n(Seq:D $s, Int:D $n --> List) {
    return [] if $n == 0;
    $s[^$n].List
}

# -------------------------
# Engel expansion
# a_k = ceil(1/r), r <- a_k*r - 1
# -------------------------

sub engel-terms($x --> Seq) {
    my $r = $x;
    gather loop {
        last if $r == 0;
        my $a = (1 / $r).ceiling;
        take $a.Int;
        $r = $a * $r - 1;
    }
}

# -------------------------
# Pierce expansion (alternating Engel)
# Greedy with floor:
# a_k = floor(1/r), r <- 1 - a_k*r
# Matches example for 1/sqrt(2): [1,3,8,33,35]
# -------------------------

sub pierce-terms($x --> Seq) {
    my $r = $x;
    gather loop {
        last if $r == 0;
        my $a = (1 / $r).floor;
        $a = 1 if $a < 1; # keep factors sensible
        take $a.Int;
        $r = 1 - $a * $r;
    }
}

# -------------------------
# Sylvester expansion (Egyptian greedy)
# a_k = ceil(1/r), r <- r - 1/a_k
# -------------------------

sub sylvester-terms($x --> Seq) {
    my $r = $x;
    gather loop {
        last if $r == 0;
        my $a = (1 / $r).ceiling;
        take $a.Int;
        $r -= 1 / $a;
    }
}

# -------------------------
# Cantor expansion (factorial number system digits)
# For n >= 0 Int:
# a_k = n mod (k+1), n <- n div (k+1), k starting at 1
# -------------------------

sub cantor-expansion($n --> List) {
    die "Cantor expansion expects a non-negative Int" unless $n ~~ Int && $n >= 0;
    my $m = $n;
    my @digits;
    my $k = 1;
    while $m > 0 {
        my $base = $k + 1;
        @digits.push: ($m % $base);
        $m div= $base;
        $k++;
    }
    @digits.List
}

# -------------------------
# Cantor product expansion (greedy)
# Choose a_k so that product P_k = Π (1 + 1/a_i) approaches x:
# Let y = x / P; choose a = floor(1/(y-1)) + 1 (for y>1).
# This matches provided examples for π.
# -------------------------

sub cantor-product-terms($x --> Seq) {
    my $target = $x;
    my $P = $target ~~ Rat|FatRat|Int ?? 1.FatRat !! 1e0;

    gather loop {
        # For exact rationals: can terminate when P == target.
        last if $P == $target;

        my $y = $target / $P;
        # If y <= 1, we cannot proceed with the standard greedy step.
        last if $y <= 1;

        my $a = (1 / ($y - 1)).floor + 1;
        $a = 1 if $a < 1;

        take $a.Int;
        $P *= (1 + 1 / $a);
    }
}

# -------------------------
# Lüroth / Lueroth expansion
# a_k = ceil(1/r)
# r <- a_k*(a_k-1)*r - (a_k-1)
#
# For exact rationals, either terminates or becomes periodic; periodic is returned as:
#   [p, [a1..ap]]
# -------------------------

sub lueroth-terms($x --> Seq:D) {
    my $r = $x;
    return do gather loop {
        last if $r == 0;
        my $a = (1 / $r).ceiling;
        take $a.Int;
        $r = $a * ($a - 1) * $r - ($a - 1);
    }
}

sub lueroth-all-rational(FatRat:D $x --> Mu) {
    my FatRat $r = $x;
    my %seen;   # remainder -> index
    my @seq;

    loop {
        return @seq.List if $r == 0;

        if %seen{$r}:exists {
            my $start = %seen{$r};
            my @period = @seq[$start .. *];
            return [ @period.elems, @period.List ];
        }

        %seen{$r} = @seq.elems;

        my $a = (1 / $r).ceiling.Int;
        @seq.push: $a;
        $r = $a * ($a - 1) * $r - ($a - 1);
    }
}

# -------------------------
# Zeckendorf representation (0/1 list over Fibonacci numbers)
# Greedy with nonconsecutive Fibonacci restriction.
# Returns bits aligned to ascending Fibonacci list starting [1,2,3,5,...]
# -------------------------

sub zeckendorf($n --> List) {
    die "Zeckendorf representation expects a non-negative Int" unless $n ~~ Int && $n >= 0;

    return [0] if $n == 0;

    my @fib = 1, 2;
    while @fib[*-1] <= $n {
        @fib.push: @fib[*-1] + @fib[*-2];
    }
    @fib.pop if @fib[*-1] > $n; # largest <= n

    my $rem = $n;
    my @bits = 0 xx @fib.elems;

    my $i = @fib.elems - 1;
    while $i >= 0 {
        if @fib[$i] <= $rem {
            @bits[$i] = 1;
            $rem -= @fib[$i];
            $i -= 2; # skip adjacent
        } else {
            $i -= 1;
        }
    }

    return @bits.reverse.List;
}

# -------------------------
# Oppenheim (placeholder in this initial module)
# Spec requires: number-expansion(x, n, r, s, p, q, "Oppenheim")
# -------------------------

multi sub number-expansion($x, Int:D $n, Int:D $r, Int:D $s, Int:D $p, Int:D $q, Str:D $t --> List) {
    my $type = normalize-type($t);
    die "This multi is only for Oppenheim expansion" unless $type eq 'Oppenheim';
    die "Oppenheim expansion is not implemented in this initial module yet (r=$r s=$s p=$p q=$q).";
}

# -------------------------
# Inverses / reconstruction
# -------------------------

# Engel:
# x ≈ Σ_{k>=1} 1 / Π_{i<=k} a_i
sub from-engel(@a --> Mu) {
    my $prod = 1.FatRat;
    my $sum  = 0.FatRat;
    for @a -> $ai {
        my $a = $ai.Int;
        die "Engel terms must be positive integers" if $a <= 0;
        $prod *= $a;
        $sum  += 1 / $prod;
    }
    return $sum;
}

# Sylvester:
# x ≈ Σ 1/a_i
sub from-sylvester(@a --> Mu) {
    my $sum = 0.FatRat;
    for @a -> $ai {
        my $a = $ai.Int;
        die "Sylvester terms must be positive integers" if $a <= 0;
        $sum += 1 / $a;
    }
    return $sum;
}

# Cantor:
# n = Σ a_k * k! , where @digits = [a1,a2,...] corresponds to 1!,2!,...
sub from-cantor(@digits --> Mu) {
    my $n = 0;
    my $fact = 1;
    my $k = 1;
    for @digits -> $d {
        $fact *= $k;               # k!
        $n += $d.Int * $fact;
        $k++;
    }
    return $n;
}

# CantorProduct:
# x ≈ Π (1 + 1/a_i)
sub from-cantor-product(@a --> Mu) {
    my $p = 1.FatRat;
    for @a -> $ai {
        my $a = $ai.Int;
        die "CantorProduct terms must be positive integers" if $a <= 0;
        $p *= (1 + 1 / $a);
    }
    return $p;
}

# Lueroth:
# x = 1/a1 + 1/(a1(a1-1)a2) + 1/(a1(a1-1)a2(a2-1)a3) + ...
# If periodic form is [p, [..]] we just reconstruct one period (useful for quick checks),
# otherwise we reconstruct from the finite prefix given.
sub from-lueroth(@terms --> Mu) {
    my @a = @terms;

    # periodic encoding: [p, [a1..ap]]
    if @a.elems == 2 && @a[0] ~~ Int && @a[1] ~~ Positional {
        @a = @a[1].List;
    }

    my $den = 1.FatRat;
    my $sum = 0.FatRat;

    for @a.kv -> $i, $ai {
        my $a = $ai.Int;
        die "Lueroth terms must be integers >= 2" if $a < 2;

        if $i == 0 {
            $den = $a;
            $sum += 1 / $den;
        } else {
            # multiply by previous (a_{i-1}-1) and current a_i
            # But easiest: keep running denominator per definition:
            # denom_k = a1(a1-1) a2(a2-1) ... a_{k-1}(a_{k-1}-1) a_k
            # We already had denom_{k-1}; update by *(a_{k-1}-1)*a_k
            my $prev = @a[$i-1].Int;
            $den *= ($prev - 1) * $a;
            $sum += 1 / $den;
        }
    }

    return $sum;
}

sub from-zeckendorf(@bits --> Mu) {
    my @pos = @bits.reverse.grep(1, :k) <<+>> 2;
    return fibonacci(@pos).sum;
}