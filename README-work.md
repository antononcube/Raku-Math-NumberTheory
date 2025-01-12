# Math::NumberTheory

Raku package with Number theory functions.

The function names and features follow the 
[Number theory functions of Wolfram Language](http://reference.wolfram.com/language/guide/NumberTheory.html).


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

### Factor integers


Gives a list of the prime factors of an integer argument, together with their exponents:

```perl6
use Math::NumberTheory;
factor-integer(factorial(20))
```

**Remark:** By default `factor-integer` uses [Pollard's Rho algorithm](https://en.wikipedia.org/wiki/Pollard's_rho_algorithm) --
specified with `method => 'rho'` --
as implemented at RosettaCode, [RC1], with similar implementations provided by "Prime::Factor", [SSp1], and "Math::Sequences", [RCp1].

Do partial factorization, pulling out at most `k` distinct factors:

```perl6
factor-integer(factorial(20), 3, method => 'trial')
```


### Chinese reminders

Data:

```perl6
my @data = 931074546, 117172357, 482333642, 199386034, 394354985;
```

Keys:

```perl6
#my @keys = random-prime(10**9 .. 10**12, @data.elems);
my @keys = 274199185649, 786765306443, 970592805341, 293623796783, 238475031661;
```

**Remark:** Using these larger keys is also a performance check.


Encrypted data:

```perl6
my $encrypted = chinese-reminder(@data, @keys);
```

Decrypted:

```perl6
my @decrypted = @keys.map($encrypted mod *);
```

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
