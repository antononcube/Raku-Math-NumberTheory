use v6.d;

unit module Math::NumberTheory::Constants;

# Using N[Sqrt[5], 100] in Wolfram Language
our constant $sqrt5 = 2.236067977499789696409173668731276235440618359611525724270897245410520925637804899414414408378782275.FatRat;

# Fibonacci 401 and 400
constant $fibonacci401 = 284812298108489611757988937681460995615380088782304890986477195645969271404032323901.FatRat;
constant $fibonacci400 = 176023680645013966468226945392411250770384383304492191886725992896575345044216019675.FatRat;

#| Golden ratio (phi)
our constant $phi = $fibonacci401 / $fibonacci400;
#our constant \ϕ is export = (1.FatRat + $sqrt5) / 2.FatRat;
sub golden-ratio(Bool:D :pre(:$pre-computed) = True) is export {
    return $pre-computed ?? $phi !! (1 + sqrt(5)) / 2;
}