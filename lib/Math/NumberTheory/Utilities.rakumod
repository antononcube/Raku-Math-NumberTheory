use v6.d;

#| Spiral square lattice
#| C<$n> -- Square side size.
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
sub triangle-matrix(Int $k where $k mod 2 == 1, :$missing-value = 0, Bool:D :d(:$dataset) = False) is export {
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