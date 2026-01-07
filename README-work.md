# Math::NumberTheory

Raku package with Number theory functions.

The function names and features follow the 
[Number theory functions of Wolfram Language](http://reference.wolfram.com/language/guide/NumberTheory.html).

**Remark:**  Raku has some nice built-in Number theory functions, like, `base`, `gcd`, `mod`, `polymod`, `expmod`, `is-prime`. 
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

## GCD

The infix operator `gcd` for calculating the Greatest Common Divisor (GCD) 
is extended to work with rational numbers and [Gaussian integers](https://en.wikipedia.org/wiki/Gaussian_integer):

### Rationals

For rational numbers `r1` and `r2`, `r1 gcd r2` gives the greatest rational number `r` for which `r1/r` and `r2/r` are integers.

```perl6
use Math::NumberTheory;
<1/3> gcd <2/5> gcd <1/7>
==> {.raku}()
```

### Gaussian integers

GCD for two Gaussian integers (complex numbers with integer real and imaginary parts):

```perl6
(10 + 15i) gcd (-3 + 2i)
```

```perl6
105 gcd (7 + 49i)
```

Here is verification of the latter:

```perl6
say 105 / (7 + 14i);
say (7 + 49i) / (7 + 14i);
```

### Prime number testing

The built-in sub `is-prime` is extended to work with [Gaussian integers](https://en.wikipedia.org/wiki/Gaussian_integer):

```perl6
say is-prime(2 + 1i);
say is-prime(5, :gaussian-integers);
```

Another extension is threading over lists of numbers:

```perl6
is-prime(^6)
```

### Factor integers

Gives a list of the prime factors of an integer argument, together with their exponents:

```perl6
factor-integer(factorial(20))
```

**Remark:** By default `factor-integer` uses [Pollard's Rho algorithm](https://en.wikipedia.org/wiki/Pollard's_rho_algorithm) --
specified with `method => 'rho'` --
as implemented at RosettaCode, [RC1], with similar implementations provided by "Prime::Factor", [SSp1], and "Math::Sequences", [RCp1].

Do partial factorization, pulling out at most `k` distinct factors:

```perl6
factor-integer(factorial(20), 3, method => 'trial')
```


### Chinese remainders

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
my $encrypted = chinese-remainder(@data, @keys);
```

Decrypted:

```perl6
my @decrypted = @keys.map($encrypted mod *);
```

### Modular exponentiation and modular inversion

The sub `power-mod` extends the built-in sub `expmod`.
The sub `modular-inverse` is based on `power-mod`.

`expmod` gives an error and no result when the 1st argument cannot be inverted with the last argument:

```perl6
expmod(30, -1, 12)
```

`power-mod` returns `Nil`:

```perl6
power-mod(30, -1, 12).defined
```


### Number-base related

There are several subs that provide functionalities related to number systems representation.

For example, here we find the digit-breakdown of $100!$:

```perl6
100.&factorial.&digit-count
```

Here is an example of using `real-digits`:

```perl6
real-digits(123.55555)
```

Non-integer bases can be also used:

```perl6
my $r = real-digits(π, ϕ);
```

Here we recover $\pi$ from the Golden ratio representation obtained above:

```perl6
$r.head.kv.map( -> $i, $d { $d * ϕ ** ($r.tail - $i - 1)  }).sum.round(10e-12);
```

The sub `phi-number-system` can be used to compute
[Phi number system](https://mathworld.wolfram.com/PhiNumberSystem.html)
representations:

```perl6
.say for (^7)».&phi-number-system
```

**Remark:** Because of the approximation of the Golden ratio used in “Math::NumberTheory”,
in order to get exact integers from phi-digits we have to round using small multiples of 10. 

-------

## CLI

The package provides the Command Line Interface (CLI) script `number-theory`. Here is its usage note:

```shell
number-theory --help
```

The script takes proper sub names as a first argument or their "conversational" form. For example, these two commands invoke the same sub:

```shell
number-theory is-happy-number 2026
```

Using ranges:

```shell
number-theory random-prime 400..440 6
```

-------

## TODO

- [ ] TODO Implementation
  - [X] DONE Gaussian integers GCD
  - [X] DONE Gaussian integers factorization
  - [X] DONE Moebius Mu function, Liouville lambda function 
    - [X] DONE Integers
    - [X] DONE Gaussian integers
  - [X] DONE Square-free test
    - [X] DONE Integers
    - [X] DONE Gaussian integers
  - [X] DONE Rational numbers GCD
  - [X] DONE Multiplicative Order
  - [X] DONE Gaussian integers LCM 
  - [X] DONE Rational numbers LCM
  - [X] DONE Carmichael lambda
  - [X] DONE Integer partitions
  - [X] DONE Number classes tests and retrievers
    - [X] DONE Abundant number
    - [X] DONE Deficient number
    - [X] DONE Perfect number
    - [X] DONE Happy number
    - [X] DONE Harshad number
  - [ ] TODO Sum of squares representation
  - [ ] TODO Figure out which memoization approach to use:
    - [ ] Via the package ["Memoize"](https://raku.land/zef:lizmat/Memoize)
    - [ ] Via `use experimental :cached` and `sub blah(...) is cached {...}` 
  - [X] DONE CLI 
- [ ] TODO Documentation
  - [X] DONE [Blog post about Number theory properties of 2026](https://raku-advent.blog/2025/12/22/day-22-numerically-2026-is-unremarkable-yet-happy/). 
  - [ ] TODO Blog post on first non-zero digit of 10_000!
  - [ ] TODO Videos
    - [X] DONE [Neat examples 1](https://www.youtube.com/watch?v=wXXWyRAAPvc)
    - [X] DONE [Neat examples 2](https://www.youtube.com/watch?v=6uCIoonlybk)
    - [X] DONE [Neat examples 3](https://www.youtube.com/watch?v=6uCIoonlybk)
      - See also the related posts [AA1, AAn4]. 
    - [ ] TODO Neat examples 4
    - [ ] TODO Neat examples 5

-------

## References

### Articles, blog posts, wiki-pages

[AA1] Anton Antonov, 
["Primitive roots generation trails"](https://mathematicaforprediction.wordpress.com/2025/04/08/primitive-roots-generation-trails/),
(2025),
[MathematicaForPrediction at WordPress](https://mathematicaforprediction.wordpress.com).

[AA2] Anton Antonov,
["Day 22 – Numerically 2026 Is Unremarkable Yet Happy"](https://raku-advent.blog/2025/12/22/day-22-numerically-2026-is-unremarkable-yet-happy/),
(2025),
[Raku Advent Calendar at WordPress](https://raku-advent.blog).

[RC1] Rosetta Code, [Prime decomposition](https://rosettacode.org/wiki/Prime_decomposition),
[Section "Pure Raku"](https://rosettacode.org/wiki/Prime_decomposition#Pure_Raku).

### Notebooks

[AAn1] Anton Antonov,
[Number theory neat examples Set 1](https://github.com/antononcube/RakuForPrediction-blog/blob/main/Presentations/Notebooks/Number-theory-neat-examples-Set-1.ipynb),
[*presentation notebook*](https://www.youtube.com/watch?v=wXXWyRAAPvc),
(2025),
[RakuForPrediction-blog at GitHub](https://github.com/antononcube/RakuForPrediction-blog).

[AAn2] Anton Antonov,
[Number theory neat examples Set 2](https://github.com/antononcube/RakuForPrediction-blog/blob/main/Presentations/Notebooks/Number-theory-neat-examples-Set-2.ipynb),
[*presentation notebook*](https://www.youtube.com/watch?v=sMwuGVvkLkU),
(2025),
[RakuForPrediction-blog at GitHub](https://github.com/antononcube/RakuForPrediction-blog).

[AAn3] Anton Antonov,
[Number theory neat examples Set 3](https://github.com/antononcube/RakuForPrediction-blog/blob/main/Presentations/Notebooks/Number-theory-neat-examples-Set-3.ipynb),
[*presentation notebook*](https://www.youtube.com/watch?v=6uCIoonlybk),
(2025),
[RakuForPrediction-blog at GitHub](https://github.com/antononcube/RakuForPrediction-blog).

[AAn4] Anton Antonov,
["Primitive roots generation trails"](https://community.wolfram.com/groups/-/m/t/3442027),
(2025),
[Wolfram Community](https://community.wolfram.com).

### Packages 

[RCp1] Raku Community,
[Math::Sequences Raku package](https://github.com/raku-community-modules/Math-Sequences),
(2016-2024),
[GitHub/raku-community-modules](https://github.com/raku-community-modules).

[SSp1] Stephen Schulze,
[Prime::Factor Raku package](https://github.com/thundergnat/Prime-Factor),
(2016-2023),
[GitHub/thundergnat](https://github.com/thundergnat).

### Videos

[AAv1] Anton Antonov,
["Number theory neat examples in Raku (Set 1)"](https://www.youtube.com/watch?v=wXXWyRAAPvc),
(2025),
[YouTube/@AAA4prediction](https://www.youtube.com/@AAA4prediction).

[AAv2] Anton Antonov,
["Number theory neat examples in Raku (Set 2)"](https://www.youtube.com/watch?v=sMwuGVvkLkU),
(2025),
[YouTube/@AAA4prediction](https://www.youtube.com/@AAA4prediction).

[AAv3] Anton Antonov,
["Number theory neat examples in Raku (Set 3)"](https://www.youtube.com/watch?v=6uCIoonlybk),
(2025),
[YouTube/@AAA4prediction](https://www.youtube.com/@AAA4prediction).