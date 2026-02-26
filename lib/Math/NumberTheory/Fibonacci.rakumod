use v6.d;

unit module Math::NumberTheory::Fibonacci;

use Math::NumberTheory::Constants;

#==========================================================
# Fibonacci
#==========================================================
# This can be cached either by using "Memoize" or by
#  use experimental :cached
# and declaration as
#  sub fibonacci(UInt $n) is cached {...}
#| Give the n-th Fibonacci number.
#| C<$n> --  Non-negative integer or a list of non-negative integers.
our proto sub fibonacci($n) is export {*}

multi sub fibonacci(@n) {
    return @n».&fibonacci.List;
}

multi sub fibonacci(UInt:D $n) {
    return 0 if $n == 0;
    return 1 if $n ≤ 2;
    my $k = $n div 2;
    if $n %% 2 {
        my $a = fibonacci($k);
        return $a * (2 * fibonacci($k + 1) - $a);
    }
    return fibonacci($k + 1) ** 2 + fibonacci($k) ** 2;
}

multi sub fibonacci(Int:D $n where $n < 0) {
    return fibonacci($n.abs) * (($n + 1) mod 2 ?? -1 !! 1);
}

multi sub fibonacci(Complex:D $n) {
    my $fn = golden-ratio() ** $n;
    return ($fn - cos($n * π) / $fn) / sqrt(5);
}

multi sub fibonacci(Rat:D $n) {
    return fibonacci($n + 0i);
}

#==========================================================
# Lucas-L
#==========================================================
#| Give the n-th Lucas L number.
#| C<$n> --  Non-negative integer or a list of non-negative integers.
our proto sub lucas-l($n) is export {*}

multi sub lucas-l(@n) {
    # Optimization
    if @n.all ~~ UInt:D {
        return fibonacci(@n <<->> 1) <<+>> fibonacci(@n <<+>> 1)
    }
    return @n».&lucas-l.List;
}

multi sub lucas-l(UInt:D $n) {
    return fibonacci($n-1) + fibonacci($n+1)
}

multi sub lucas-l(Int:D $n where $n < 0) {
    return lucas-l($n.abs) * ($n mod 2 ?? -1 !! 1);
}
multi sub lucas-l(Complex:D $n) {
    my $fn = golden-ratio() ** $n;
    return $fn + cos($n * π) / $fn;
}

multi sub lucas-l(Rat:D $n) {
    return lucas-l($n + 0i);
}