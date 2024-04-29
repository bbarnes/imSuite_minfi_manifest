#! /usr/bin/perl -w

use Data::Dumper;
 
package Foo;
sub new {bless {'a' => 1, 'b' => sub { return "foo" }}, $_[0]};
 
package Fuz;                       # a weird REF-REF-SCALAR object
sub new {bless \($_ = \ 'fu\'z'), $_[0]};
 
package main;
$foo = Foo->new;
$fuz = Fuz->new;
$boo = [ 1, [], "abcd", \*foo,
{1 => 'a', 023 => 'b', 0x45 => 'c'},
    \\"p\q\'r", $foo, $fuz];
 
########
# simple usage
########
 
$bar = eval(Dumper($boo));
print($@) if $@;
print Dumper($boo), Dumper($bar);  # pretty print (no array indices)
 
$Data::Dumper::Terse = 1;        # don't output names where feasible
$Data::Dumper::Indent = 0;       # turn off all pretty print
print Dumper($boo), "\n";
 
$Data::Dumper::Indent = 1;       # mild pretty print
print Dumper($boo);
 
$Data::Dumper::Indent = 3;       # pretty print with array indices
print Dumper($boo);
 
$Data::Dumper::Useqq = 1;        # print strings in double quotes
print Dumper($boo);
 
$Data::Dumper::Pair = " : ";     # specify hash key/value separator
print Dumper($boo);
 
 
########
# recursive structures
########
 
@c = ('c');
$c = \@c;
$b = {};
$a = [1, $b, $c];
$b->{a} = $a;
$b->{b} = $a->[1];
$b->{c} = $a->[2];
print Data::Dumper->Dump([$a,$b,$c], [qw(a b c)]);
 
 
$Data::Dumper::Purity = 1;         # fill in the holes for eval
print Data::Dumper->Dump([$a, $b], [qw(*a b)]); # print as @a
print Data::Dumper->Dump([$b, $a], [qw(*b a)]); # print as %b
 
 
$Data::Dumper::Deepcopy = 1;       # avoid cross-refs
print Data::Dumper->Dump([$b, $a], [qw(*b a)]);
 
 
$Data::Dumper::Purity = 0;         # avoid cross-refs
print Data::Dumper->Dump([$b, $a], [qw(*b a)]);
 
########
# deep structures
########
 
$a = "pearl";
$b = [ $a ];
$c = { 'b' => $b };
$d = [ $c ];
$e = { 'd' => $d };
$f = { 'e' => $e };
print Data::Dumper->Dump([$f], [qw(f)]);
 
$Data::Dumper::Maxdepth = 3;       # no deeper than 3 refs down
print Data::Dumper->Dump([$f], [qw(f)]);
 
 
########
# object-oriented usage
########
 
$d = Data::Dumper->new([$a,$b], [qw(a b)]);
$d->Seen({'*c' => $c});            # stash a ref without printing it
$d->Indent(3);
print $d->Dump;
$d->Reset->Purity(0);              # empty the seen cache
print join "----\n", $d->Dump;
 
 
########
# persistence
########
 
package Foo;
sub new { bless { state => 'awake' }, shift }
sub Freeze {
    my $s = shift;
    print STDERR "preparing to sleep\n";
    $s->{state} = 'asleep';
    return bless $s, 'Foo::ZZZ';
}
 
package Foo::ZZZ;
sub Thaw {
    my $s = shift;
    print STDERR "waking up\n";
    $s->{state} = 'awake';
    return bless $s, 'Foo';
}
 
package main;
use Data::Dumper;
$a = Foo->new;
$b = Data::Dumper->new([$a], ['c']);
$b->Freezer('Freeze');
$b->Toaster('Thaw');
$c = $b->Dump;
print $c;
$d = eval $c;
print Data::Dumper->Dump([$d], ['d']);
 
 
########
# symbol substitution (useful for recreating CODE refs)
########
 
sub foo { print "foo speaking\n" }
*other = \&foo;
$bar = [ \&other ];
$d = Data::Dumper->new([\&other,$bar],['*other','bar']);
$d->Seen({ '*foo' => \&foo });
print $d->Dump;
 
 
########
# sorting and filtering hash keys
########
 
$Data::Dumper::Sortkeys = \&my_filter;
my $foo = { map { (ord, "$_$_$_") } 'I'..'Q' };
my $bar = { %$foo };
my $baz = { reverse %$foo };
print Dumper [ $foo, $bar, $baz ];
 
sub my_filter {
    my ($hash) = @_;
    # return an array ref containing the hash keys to dump
    # in the order that you want them to be dumped
    return [
      # Sort the keys of %$foo in reverse numeric order
        $hash eq $foo ? (sort {$b <=> $a} keys %$hash) :
      # Only dump the odd number keys of %$bar
        $hash eq $bar ? (grep {$_ % 2} keys %$hash) :
      # Sort keys in default order for all other hashes
        (sort keys %$hash)
	];
}

## End of file
