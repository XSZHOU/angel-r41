#! /usr/bin/perl -w

$line=<>;
$line =~ s/\t/\ \ \ \ /g;
$line =~ s/\ *$//;
$lastline = $line;

while ($line=<>) {
    $line =~ s/\t/\ \ \ \ /g;
    $line =~ s/\ *$//;

    print $lastline;

    $lastline = $line;
}

unless ($lastline =~ /^\n$/) {
    if ($lastline =~ /\n$/) {
        print $lastline;
    } else {
        print "$lastline\n";
    }
}
