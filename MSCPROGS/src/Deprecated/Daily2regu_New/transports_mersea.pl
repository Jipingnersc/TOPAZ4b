#!/usr/bin/perl
use lib qw(/home/parallab/nersc/knutali/Perl/lib,
           /home/parallab/nersc/knutali/Perl/lib/site_perl/);
use warnings;
use strict 'vars';
use NetCDF;
use IO::File;

my $bull_date;
my $bull_date2;
my $forc_date;
my ($dimlen, $dimname);
my ($ncid, $ncfile);
my ($ncid_trvar,$ncid_tridvar);
my ($trvar,$tridvar);
my ($dimid,$attid);
my ($status);
my $id;
my $filehandle;
my $fname;
my $line;
my $trval;
my $count;

# Pass bulletin and forecast as arguments
$bull_date=$ARGV[0];
$forc_date=$ARGV[1];

$bull_date2=$bull_date;
$bull_date=~s/f//;


# stuff
$dimlen=136; # Number of transports
$dimname="transportid";
$trvar="transport";
$tridvar="transportid";


# Create file name
$ncfile="topaz_V1_mersea_nat_transport_b${bull_date}_f${forc_date}9999.nc";

# Open Netcdf file
$ncid=NetCDF::create($ncfile, NetCDF::CLOBBER);

# Global attributes
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "title" , NetCDF::CHAR,
   "TOPAZ Daily average");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "institution" , NetCDF::CHAR, 
   "NERSC, Edvard Griegs Vei 3A, Bergen, Norway");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "history" , NetCDF::CHAR, 
   "Originally Created from TOPAZ daily averages for MERSEA") ;
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "references" , NetCDF::CHAR, 
   "http://www.mersea.eu.org ; http://topaz.nersc.no" );
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "bulletin_type" , NetCDF::CHAR, 
   "Hindcast");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "bulletin_date" , NetCDF::CHAR, 
   "${bull_date}");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "field_date" , NetCDF::CHAR, 
   "${forc_date}");



# Create dimension
$dimid = NetCDF::dimdef($ncid, $dimname, $dimlen );


# Define variables and attributes for
# 1) Transport 
$ncid_trvar= NetCDF::vardef($ncid, $trvar, NetCDF::FLOAT, $dimid);
$attid = NetCDF::attput($ncid, $ncid_trvar, "long_name", NetCDF::CHAR, "section transport");
$attid = NetCDF::attput($ncid, $ncid_trvar, "axis", NetCDF::CHAR, "X");
$attid = NetCDF::attput($ncid, $ncid_trvar, "standard_name", NetCDF::CHAR, "transport");
$attid = NetCDF::attput($ncid, $ncid_trvar, "units", NetCDF::CHAR, "Sverdrup");

# 2) Transport id
$ncid_tridvar= NetCDF::vardef($ncid, $tridvar, NetCDF::SHORT, $dimid);
$attid = NetCDF::attput($ncid, $ncid_trvar, "long_name", NetCDF::CHAR, "transport Identification");


# End definition
$status = NetCDF::endef($ncid);



#Put in ze data -- this assumes all files are of the form "TRP001.dat"
$count=0;
while (<*TRP*.dat>){

   # Transport ID is in the file name -- This retrieves the relevant part for 
   # TOPAZ file names
   $fname=$_;
   $id=$_;
   $id=~s/.*TRP//;
   $id=~s/_b.*//;

   # TOPAZ files have three lines. the last line contains the actual transport
   $filehandle = new IO::File;
   $filehandle->open("<$fname") or die "Could not open file $fname";

   $filehandle->getline();
   $filehandle->getline();
   $line=$filehandle->getline();
   $trval=((split " ","$line")[1]);
   $filehandle->close;

   print "$trval\n";

   $status = NetCDF::varput($ncid, $ncid_tridvar, $count, 1, $id   );
   $status = NetCDF::varput($ncid, $ncid_trvar  , $count, 1, $trval);

   $count=$count+1;
}

NetCDF::close($ncid) == 0 || die "Couldn't close netCDF file\n";

