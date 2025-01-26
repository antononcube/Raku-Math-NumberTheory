# Math::NumberTheory

Raku package with Number theory functions.

The function names and features follow the 
[Number theory functions of Wolfram Language](http://reference.wolfram.com/language/guide/NumberTheory.html).

**Remark:**  Raku has some nice built-in Number theory functions, like, `base`, `mod`, `polymod`, `expmod`, `is-prime`. 
They somewhat lack generality, hence their functionalities are extended with this package. 
For example, `is-prime` works with lists and [Gaussian integers](https://en.wikipedia.org/wiki/Gaussian_integer).

-------

## Installation

From [Zef ecosystem](https://raku.land):

```
zef install Math::NumberTheory
```

From GitHub:

```
zef install https://github.com/antononcube/Raku-Math-NumberTheory
```

-------

## Usage examples

### Prime number testing

The built-in sub `is-prime` is extended to work with [Gaussian integers](https://en.wikipedia.org/wiki/Gaussian_integer):

```perl6
use Math::NumberTheory;
say is-prime(2 + 1i);
say is-prime(5, :gaussian-integers);
```
```
# True
# False
```

Another extension is threading over lists of numbers:

```perl6
is-prime(^6)
```
```
# (False False True True False True)
```

### Factor integers

Gives a list of the prime factors of an integer argument, together with their exponents:

```perl6
factor-integer(factorial(20))
```
```
# [(2 18) (3 8) (5 4) (7 2) (11 1) (13 1) (17 1) (19 1)]
```

**Remark:** By default `factor-integer` uses [Pollard's Rho algorithm](https://en.wikipedia.org/wiki/Pollard's_rho_algorithm) --
specified with `method => 'rho'` --
as implemented at RosettaCode, [RC1], with similar implementations provided by "Prime::Factor", [SSp1], and "Math::Sequences", [RCp1].

Do partial factorization, pulling out at most `k` distinct factors:

```perl6
factor-integer(factorial(20), 3, method => 'trial')
```
```
# [(2 18) (3 8) (5 4)]
```


### Chinese reminders

Data:

```perl6
my @data = 931074546, 117172357, 482333642, 199386034, 394354985;
```
```
# [931074546 117172357 482333642 199386034 394354985]
```

Keys:

```perl6
#my @keys = random-prime(10**9 .. 10**12, @data.elems);
my @keys = 274199185649, 786765306443, 970592805341, 293623796783, 238475031661;
```
```
# [274199185649 786765306443 970592805341 293623796783 238475031661]
```

**Remark:** Using these larger keys is also a performance check.


Encrypted data:

```perl6
my $encrypted = chinese-reminder(@data, @keys);
```
```
# 6681669841357504673192908619871066558177944924838942629020
```

Decrypted:

```perl6
my @decrypted = @keys.map($encrypted mod *);
```
```
# [931074546 117172357 482333642 199386034 394354985]
```

### Modular exponentiation and modular inversion

The sub `power-mod` extends the built-in sub `expmod`.
The sub `modular-inverse` is based on `power-mod`.

`expmod` gives an error and no result when the 1st argument cannot be inverted with the last argument:

```perl6
expmod(30, -1, 12)
```
```
#ERROR: Error in mp_exptmod: Value out of range
# Nil
```

`power-mod` returns `Nil`:

```perl6
power-mod(30, -1, 12).defined
```
```
# False
```


### Number-base related

There are several subs that provide functionalities related to number systems representation.

For example, here we find the digit-breakdown of $100!$:

```perl6
100.&factorial.&digit-count
```
```
# {0 => 30, 1 => 15, 2 => 19, 3 => 10, 4 => 10, 5 => 14, 6 => 19, 7 => 7, 8 => 14, 9 => 20}
```

Here is an example of using `real-digits`:

```perl6
real-digits(123.55555)
```
```
# ([1 2 3 5 5 5 5 5] 3)
```

Non-integer bases can be also used:

```perl6
my $r = real-digits(π, ϕ);
```
```
# ([1 0 0 0 1 0 0 1 0 1 0 1 0 0 1 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 1 0 1 0 1 0 1] 3)
```

Here we recover $\pi$ from the Golden ratio representation obtained above:

```perl6
$r.head.kv.map( -> $i, $d { $d * ϕ ** ($r.tail - $i - 1)  }).sum.round(10e-12);
```
```
# 3.1415926535899996
```

The sub `phi-number-system` can be used to compute
[Phi number system](https://mathworld.wolfram.com/PhiNumberSystem.html)
representations:

```perl6
.say for (^7)».&phi-number-system
```
```
# ()
# (0)
# (1 -2)
# (2 -2)
# (2 0 -2)
# (3 -1 -4)
# (3 1 -4)
```

**Remark:** Because of the approximation of the Golden ratio used in “Math::NumberTheory”,
in order to get exact integers from phi-digits we have to round using small multiples of 10. 

-------

## TODO

- [ ] TODO Implementation
  - [ ] TODO Gaussian integers factorization
  - [ ] TODO Integer partitions
  - [ ] TODO Square-free
  - [ ] TODO Sum of squares representation
  - [ ] TODO CLI
- [ ] TODO Documentation
  - [ ] TODO Blog post on first non-zero digit of 10_000!
  - [ ] TODO Videos
    - [ ] TODO Neat examples 1
    - [ ] TODO Neat examples 2
    - [ ] TODO Neat examples 3

-------

## References

[RC1] Rosetta Code, [Prime decomposition](https://rosettacode.org/wiki/Prime_decomposition),
[Section "Pure Raku"](https://rosettacode.org/wiki/Prime_decomposition#Pure_Raku).

[RCp1] Raku Community,
[Math::Sequences Raku package](https://github.com/raku-community-modules/Math-Sequences),
(2016-2024),
[GitHub/raku-community-modules](https://github.com/raku-community-modules).

[SSp1] Stephen Schulze,
[Prime::Factor Raku package](https://github.com/thundergnat/Prime-Factor),
(2016-2023),
[GitHub/thundergnat](https://github.com/thundergnat).
