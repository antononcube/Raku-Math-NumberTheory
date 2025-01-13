use v6.d;

#| Spiral square lattice
#| C<$n> -- Square side size.
sub spiral-lattice(UInt:D $n) is export {
    my @matrix = (0 xx $n).Array xx $n ;
    my ($row, $col, $num) = (0, 0, 1);
    my @directions = ([0, 1], [1, 0], [0, -1], [-1, 0]);
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
    return $n² + 1 <<->> @matrix;
}