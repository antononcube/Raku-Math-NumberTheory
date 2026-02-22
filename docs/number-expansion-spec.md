# Number expansions

The sub `number-expansion` can express a number as an expansion using any of several methods.

### Documentation

#### Usage

- `number-expansion(x, n, t)`
  - generates a list of the first n terms in the series representation of x **for the chosen expansion type *t*.

- `number-expansion(x, n, t)`
  - generates a list of all the terms that can be obtained using arbitrary-precision arithmetic.

#### Details & Options

- "Lüroth" can also be written as "Lueroth" when passed as an argument.

- The Engel expansion representation $\{a_1, a_2, ...\}$ corresponds to the expression $1/a_1 + 1/(a_1 a_2) + 1/(a_1 a_2 a_3)+...$

- The Pierce expansion (alternating Engel expansion) representation $\{a_1, a_2, ...\}$ corresponds to the expression $1/a_1 - 1/(a_1 a_2) + 1/(a_1 a_2 a_3)+...$

- The Sylvester expansion representation $\{a_1, a_2, ...\}$ corresponds to the expression $1/a_1 + 1/a_2 + ...$.

- The Cantor expansion representation $\{a_1, a_2, ...\}$ corresponds to the expression $a_1 * 1! + a_2 * 2! + a_3 * 3!...$.

- The Cantor product expansion $\{a_1, a_2, ...\}$ representation corresponds to the expression $(1+1/a_1)(1+1/a_2)(1+1/a_3)...$;

- The Lüroth expansion representation $\{a_1, a_2, ...\}$ corresponds to the expression $1/a_1 + 1/(a_1 (a_1 - 1) a_2) +1/(a_1 (a_1 - 1) a_2 (a_2 - 1) a_3) + ...$.

- The Oppenheim expansion representation $\{d_1, d_2, ...\}$ corresponds to the expression $x = \frac{1}{d_1} + \sum_{n=1}^{\infty} \left( \prod_{i=1}^{n} \frac{a_i(d_i)}{b_i(d_i)} \right) \frac{1}{d_{n+1}}$.

- The Oppenheim expansion requires explicit specification of the constants a and b as defined in the original paper by A. Oppenheim.
  `number-expansion(x, n, r, s, p, q, "Oppenheim")` generates a list of the first $n$ terms in the Oppenheim series representation of $x$ where $a_i = r + s * d_i$ and $b_i = p + q + d_i, i=1,2,...,n$.

- The Zeckendorf representation gives the $0-1$ list that indicates the unique nonconsecutive Fibonacci numbers that sum to the non-negative integer $n$.

- The x can be either an exact or an inexact number.

- For exact numbers, `number-expansion(x, t)` can be used if x is rational.

- Since irrational numbers always yield an infinite sequence, the number of terms has to be specified explicitly.

- Since the series expansion representation for a rational number has only a limited number of terms, `number-expansion(x, n, t)` may yield a list with fewer than n elements in this case.

- Lüroth expansion always gives a terminating sequence, or an infinite periodic sequence for rational numbers. The latter is represented as $\{p, \{a_1, a_2, ..., a_p\}\}$, where $p$ is the periodicity.

- The function `from-number-expansion` reconstructs a number from the result of `number-expansion`.

### Examples

#### Basic Examples

10 terms in the Engel expansion of π:

```raku
number-expansion(Pi, 10, "Engel")
```

```
[1, 1, 1, 8, 8, 17, 19, 300, 1991, 2492]
```

#### Scope

Expand a rational number:

```raku
number-expansion(11/18, "Lueroth")
```

```
[2, 5, 3, 2, 3]
```

Compute the original value using the definition of the expansion:

```raku
1/2 + 1/(2*1*5) + 1/(2*1*5*4*3) + 1/(2*1*5*4*3*2*2) + 1/(2*1*5*4*3*2*2*1*3)
```

```
11/18
```

Engel expansion of the number $1.175$:

```raku
number-expansion(1.175, "Engel")
```

```
[1, 6, 20]
```


First 5 terms of the Pierce expansion of the number $1/\sqrt 2$:

```raku
number-expansion(1/sqrt(2), 5, "Pierce")
```

```
[1, 3, 8, 33, 35]
```


Sylvester expansion of the rational number $3/19$:

```raku
number-expansion(3/19, "Sylvester")
```

```
[7, 67]
```


Cantor expansion of the number $384$:

```raku
number-expansion(384, "Cantor")
```

```
[0, 0, 0, 1, 3]
```

First 5 terms of the Cantor product expansion of the irrational number $\pi$:

```raku
number-expansion(pi, 5, "CantorProduct")
```

```
[1, 2, 22, 600, 1800856, 15150670259531]
```


Lüroth expansion of the rational number 5/13. It returns an infinite expansion series with periodicity $3$:

```raku
number-expansion(5/13, "Lueroth")
```

```
[3, [3, 4, 2]]
```

Lüroth expansion of the golden ration reciprocal:

```raku
number-expansion(1/golden-ratio, 10, "Lueroth")
```

```
[2, 5, 2, 3, 2, 4, 2, 2, 162, 2]
```

Oppenheim expansion of 1/π. See Details and Options for the extra parameters in the Oppenheim expansion:

```raku
number-expansion(1/pi, 10, 1, 0, 0, 1, "Oppenheim")
```

```
[4, 4, 11, 45, 70, 1111, 4423, 5478, 49340, 94388, 200677]
```

#### Applications

Fractions of the form $a\left/3^k\right.$ are conjectured to always have a finite Lüroth expansion. This can be investigated with the following:

```raku
number-expansion(7/27, "Lueroth")
```

```
(4, 9)
```

```raku
number-expansion(13/81, "Lueroth")
```

```
(7, 2, 3, 2, 2, 2, 9)
```

Whether all terms of the Sylvester series (Sylvester expansion of the number 1) are square-free is an open problem. All known terms till now are square-free. This urges us to further investigate which terms might have all square-free terms. The Sylvester expansions for `1/Pi` and `1/GoldenRatio` have square terms:

```raku
number-expansion(1/pi, 5, "Sylvester")
```

```
(4, 15, 609, 845029, 1010073215739)
```

```raku
number-expansion(1/golden-ratio, 5, "Sylvester")
```

```
(2, 9, 145, 37986, 2345721887)
```

From this, it might seem numbers greater than or equal to 1 might have all square-free terms. But the Sylvester expansion of ϕ seems to have a square term:

```raku
number-expansion[GoldenRatio, 5, "Sylvester"]
```

```
(1, 2, 9, 145, 37986)
```

#### Properties and Relations

The resource function [Fromnumber-expansion](https://resources.wolframcloud.com/FunctionRepository/resources/Fromnumber-expansion) is effectively the inverse of `number-expansion`:

```raku
number-expansion(pi, 10, "Engel")
```

```
(1, 1, 1, 8, 8, 17, 19, 300, 1991, 2492)
```

```raku
from-number-expansion(_, "Engel")
```

```
3.14159
```

#### Possible Issues

Expanding an irrational number requires specifying a finite length:

```raku
number-expansion(sqrt(2), "CantorProduct")
```

```
# error
```

```raku
number-expansion(sqrt(2), 5, "CantorProduct")
```

```
(3, 17, 577, 665857, 886731088897, 1572584048032918633353217)
```

```raku
number-expansion(pi, "Engel")
```

```
# error
```

Since the function uses arbitrary-precision arithmetic, computing a large number of terms might be very slow depending on the user's hardware.

### Source & Additional Information

#### Keywords

- Series expansion

- Real numbers

- Generalized Number Expansion

- Engel Expansion

- Pierce Expansion

- Alternating Engel Expansion

- Sylvester Expansion

- Cantor Expansion

- Cantor Product

- Luroth expansion

- Openheim expansion


#### Links

- [Wikipedia--Engel expansion](https://en.wikipedia.org/wiki/Engel_expansion)

- [Knopfmacher Expansions in Number Theory](https://www.tandfonline.com/doi/pdf/10.1080/16073606.2001.9639227)

- [Normal Numbers With Respect to the Cantor Series Expansion](http://rave.ohiolink.edu/etdc/view?acc_num=osu1274431587)

- [Approximation of Real Numbers by Rationals via Series Expansions](http://ac.inf.elte.hu/Vol_018_1999/089.pdf)

- [Edifice of the Real Numbers by Alternating Series](http://www.ijma.info/index.php/ijma/article/view/1587)

- [`NumberExpansion` at Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/resources/NumberExpansion)
