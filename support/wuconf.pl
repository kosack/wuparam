#!/usr/bin/perl
# wuconf.pl 
#
# example perl script for extracting information from wuparam
# configuration files.


while (<STDIN>) {

    chomp;      # remove the end-of-line character
    s/\#.*//;   # strip off comments

    # split the line at the '=' sign and skip any lines that fail
    # (ones that don't have key=value format)

    ($key,$value) = split("=") or next ;    

    # now do stuff with the key and value for example:

    print "KEY: '$key' VALUE: '$value' \n" ;

    # example of pattern matching: /pattern/i, means look for the
    # string 'pattern', ignoring case.  Using a regexp is safer than
    # comparing the strings since people may have whitespace
    # characters or other junk in the file. Of course, /pr/i will
    # match any key that has the letters "pr" in it, so you may want
    # to be more specific.

    if ( $key =~ /^pr/i ) {
	
	($on,$off,$n2,$utdate) = split(" ",$value);

	print "RUN PAIR: [$on,$off] N2:$n2 UT:$utdate\n";

    }
    

}
