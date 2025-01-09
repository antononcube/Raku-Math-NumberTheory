#!/usr/bin/env raku
use v6.d;

use Test;

use lib <. lib>;
use Math::NumberTheory;

## 1
is-deeply factor-integers(120), ((2, 3), (3, 1), (5, 1));