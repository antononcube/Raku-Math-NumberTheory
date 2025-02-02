use v6.d;

#| Spiral square lattice
#| C<$n> -- Square side size.
#| C<:l(:end-corner(:$last-at))> -- Corner to finish the spiral at.
#| C<:d(:$dataset)> -- Should a dataset be returned or not.
sub spiral-lattice(UInt:D $n, :l(:end-corner(:$last-at)) is copy = Whatever, Bool:D :d(:$dataset) = False) is export {

    if $last-at.isa(Whatever) { $last-at = 'bottom-right' }
    my @corners = <top-left top-right bottom-right bottom-left>;
    die "The value of \$last-at is expected to be Whatever or one of the strings: ⎡{@corners.join(' ')}⎦."
    unless $last-at ~~ Str:D && $last-at ∈ @corners;

    my ($row, $col, $num, @directions);
    given $last-at {
        when 'top-left' {
            ($row, $col, $num) = (0, 0, 1);
            @directions = [0, 1], [1, 0], [0, -1], [-1, 0];
        }
        when 'bottom-right' {
            ($row, $col, $num) = ($n - 1, $n - 1, 1);
            @directions = [0, -1], [-1, 0], [0, 1], [1, 0];
        }
        when 'top-right' {
            ($row, $col, $num) = (0, $n - 1, 1);
            @directions = [0, -1], [1, 0], [0, 1], [-1, 0];
        }
        when 'bottom-left' {
            ($row, $col, $num) = ($n - 1, 0, 1);
            @directions = [0, 1], [-1, 0], [0, -1], [1, 0];
        }
    }

    my @matrix = (0 xx $n).Array xx $n;
    my $dir = 0;

    while $num <= $n * $n {
        @matrix[$row][$col] = $num++;
        my $next-row = $row + @directions[$dir][0];
        my $next-col = $col + @directions[$dir][1];
        if $next-row < 0 || $next-row >= $n || $next-col < 0 || $next-col >= $n || @matrix[$next-row][$next-col] {
            $dir = ($dir + 1) % 4;
        }
        $row += @directions[$dir][0];
        $col += @directions[$dir][1];
    }

    @matrix = $n² + 1 <<->> @matrix;

    return do if $dataset {
        @matrix.map(*.kv.Hash).Array
    } else {
        @matrix
    }
}

#| Gives a triangle within a matrix.
##| C<$k> -- Number of rows of the embedding matrix. (Odd numbers only.)
##| C<:na(:$missing-value))> -- Missing value.
##| C<:d(:$dataset)> -- Should a dataset be returned or not.
sub triangle-matrix-embedding(Int:D $k, :na(:$missing-value) = 0, Bool:D :d(:$dataset) = False) is export {

    die "The first argument is expected to be a positve integer."
    unless $k > 0;

    my $ncols = (2 * $k - 1);
    my @matrix = [$missing-value xx $ncols] xx $k;
    my $start = 1;
    my $mid = $ncols div 2;

    for ^$k -> $row {
        my $offset = $row * 2;
        my $num_elements = $row * 2 + 1;

        for ^$num_elements -> $col {
            @matrix[$row][$mid - $row + $col] = $start++;
        }
    }

    return do if $dataset {
        @matrix.map(*.kv.Hash).Array
    } else {
        @matrix
    }
}

constant $golden-ratio = (1 + sqrt(5)) / 2;

constant $sunMsg = "The first argument is expected to be a positve integer or list of integers.";

#| Sunflower embedding of integers from 1 to given upper limit.
#| C<$n> -- A positive integer or a list of integers.
#| C<:&with> -- Function to apply to each integer and add the result to the corresponding record.
#| If C<WhateverCode> then no such application and addition is done.
proto sub sunflower-embedding($n, :&with = WhateverCode, Bool:D :d(:$dataset) = False) is export {*}

multi sub sunflower-embedding(Int:D $n, :&with = WhateverCode, Bool:D :d(:$dataset) = False) {
    return sunflower-embedding((1...$n).Array, :&with, :$dataset);
}
multi sub sunflower-embedding(@ints, :&with = WhateverCode, Bool:D :d(:$dataset) = False) {
    die $sunMsg unless @ints.all ~~ Int:D;

    my @sunflower = @ints.map({
        my $a = $_;
        my $angle = $a * 2 * π / $golden-ratio ** 2;
        my %res = x => sqrt($a) * cos($angle), y => sqrt($a) * sin($angle);
        if &with ~~ Callable && !&with.isa(WhateverCode) {
            %res<group> = &with($a)
        }
        %res
    });

    return do if $dataset {
        @sunflower
    } else {
        my @keys = &with ~~ Callable && !&with.isa(WhateverCode) ?? <x y group> !! <x y>;
        @sunflower.map(*{@keys})
    }
}