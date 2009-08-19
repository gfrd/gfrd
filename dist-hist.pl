#!/usr/bin/perl

use POSIX;

$pi = 3.1415;

$filename = "/storage1/bossen/devel/epdp-0.2/GF2DRad_test.distdata";
$sigma = 1e-7;
$a = 5e-7;
$num_r_bins = 10;
$num_theta_bins =20;

$r_bin_size = ($a - $sigma)/$num_r_bins;
$theta_bin_size = 2*$pi/$num_theta_bins;
$count = 0;					# counts the total number of datapoints processed

# intitializing the array
for $r (1 .. $num_r_bins) {
	for $theta (1 .. $num_theta_bins) {
		$histogram [$r][$theta] = 0;
	}
}

open (FILE, $filename) || die ("cannot open the file ");

# processing the data into bins
while (<FILE>) {
	chomp ;	# remove the last enter?
	@getallen = split (/ /);
	$r = $getallen [6];		# get the r and theta from the data
	$theta = $getallen [7];

	$r_bin_num = ceil(($r - $sigma)/$r_bin_size);
	$theta_bin_num = ceil($theta/$theta_bin_size);
	$histogram[$r_bin_num][$theta_bin_num] += 1;
#	print "$r_bin_num $theta_bin_num\n";

	$count += 1;		# increase the counter
}
close FILE;

#printing the array
for $r (1 .. $num_r_bins) {
        for $theta (1 .. $num_theta_bins) {
		$r_pos = ($r-0.5)*$r_bin_size + $sigma;
		$theta_pos = ($theta-0.5)*$theta_bin_size;
		$hight = $histogram[$r][$theta]/$count;
                print "$r_pos $theta_pos $hight\n";
        }
}


