#!/usr/bin/perl 
while(<>) {
	next if /^#.*/;
	print $_;
}
