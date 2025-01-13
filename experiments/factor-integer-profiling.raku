#!/usr/bin/env raku
use v6.d;

use Math::NumberTheory;
#use Math::Sequences::Integer :support;

my &test-func = &factor-integer;
#my &test-func = &factors;

#((2, 3), (3, 1), (5, 1));
say &test-func(120);

#say &test-func(factorial(20), 4);
{
    use Math::NumberTheory;
    say "trial-divisors 20! : ", factor-integer(factorial(20), method => 'trial');
    say "pollard-rho    20! : ", factor-integer(factorial(20), method => 'rho');
    say "pollard-rho 15227063669158801 : ", factor-integer(15227063669158801, method => 'rho');
}

my $method = 'trial';
my $n = factorial(20);
my $tstart = now;
my $res = &test-func($n);
my $res2 = $res.raku;
my $tend = now;
say "time: { $tend - $tstart } for 20! \nresult: $res2";


# Large integers
$n = 10 ** 50 + 3;
my $tstart2 = now;
$res = &test-func($n);
$res2 = $res.raku;
my $tend2 = now;
say "time: { $tend2 - $tstart2 } for 10 ** 50 + 3\nresult: $res2";