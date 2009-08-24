#!/usr/bin/perl

use POSIX;

$num_bins = 16;

$bin_size = (16)/$num_bins;
$count = 0;					# counts the total number of datapoints processed

# intitializing the array
for $r (1 .. $num_bins) {
	$histogram [$r] = 0;
}

# processing the data into bins
while (<>) {
	chomp ;	# remove the last enter?
	$r = $_;		# get the r and theta from the data
	if ($r > 0.0) {
		$log_r = -log($r)/log(10);
		$bin_num = ceil($log_r/$bin_size);
#		print "$r $log $bin_num\n";
		$histogram[$bin_num] += 1;

		$count += 1;		# increase the counter
	}
}

#printing the array
print "total processed datapoints: $count \n";
for $r (1 .. $num_bins) {
	$hight = $histogram[$r];
        print "$r $hight\n";
}


