#!/usr/bin/perl 
# animateplot - plot a bunch of images from wuparam into an animation

if (scalar(@ARGV != 1)) {
    print "USAGE: animateplot.pl <wuplot template file>\n";
    print "\twhere <wuplot template file> is a wuplot file with\n";
    print "\tthe following template variables:\n";
    print "\t\t [EXCESSIMAGE] - name of excess image\n";
    print "\t\t [SIGNIFIMAGE] - name of significance image\n\n";
    exit(1);
}

$plotfile = shift();

open(PLOTFILE, $plotfile) or die "Couldn't open '$plotfile'\n";
$/=undef;
$plottext = <PLOTFILE>;
close(PLOTFILE);
if (! (-e "Animation")) {
    mkdir("Animation") or die "Couldn't create Animation directory\n";
}
$i = 0;

foreach $file (<gt*/excess.im2d>) {
    
    $excess = $file;
    $signif = $file;
    $signif =~ s/excess/signif/;
    $outfile = sprintf("Animation/%08d",$i);

    print "Processing $file - $excess - $signif\n";

    $newplottext = $plottext;
    $newplottext =~ s/\[EXCESSIMAGE\]/$excess/g;
    $newplottext =~ s/\[SIGNIFIMAGE\]/$signif/g;

    open(OUTFILE, ">$outfile.wuplot");
    print OUTFILE $newplottext;
    close(OUTFILE);

    `wuplot -o $outfile.png -T png $outfile.wuplot`;

    $i++;

}
