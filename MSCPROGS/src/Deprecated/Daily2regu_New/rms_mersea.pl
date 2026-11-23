#!/usr/bin/perl -w
use lib qw(/home/parallab/nersc/knutali/Perl/lib,
           /home/parallab/nersc/knutali/Perl/lib/site_perl/);
use warnings;
use strict 'vars';
use NetCDF;
use IO::File;
use Date::Calc qw(Delta_Days);
use POSIX qw(ceil floor);

my $bull_date;
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
my ($region, $creg);
my ($intf, $intf_dim, $cint);
my ($RegionId, @RMSDims3D, @RMSDims2D, @RMSDims1D, $RMSSSTId, $RMSSSTName);
my ($RMSTempName, $RMSSalnName, $RMSTideName);
my ($RMSTempId  , $RMSSalnId  , $RMSTideId  );
my (@obstypes, @VarIds, @infiles,$obstype);
my ($FLDVarId, $RMSTideVarId, @Var3D, $obsind);
my ($ifile, @flds, $nrobs, $floattime, $rmsval);
my (%DimLen,%DimName, %DimId,%VarName,%VarId);
my (%FLDVarName,%FLDVarId, %obsDim, @start, @count);
my ($forc_dim, $forc_mnd, $forc_yr, $bull_dim, $bull_mnd,
   $bull_yr,$delta,@RMSvals);
my ($maxbin, @RMSnrobs, $obsindx, $bin1, $bin2,%RmsTypeVarName,%FLDVarUnit);
my ($undef, @firstregion, @firstinterval,@textdims, $subline, $verbose, %FileVarName, %Average);

$undef=-999;
#$verbose=0; # Set to 1 for verbose output
$verbose=1; # Set to 1 for verbose output

# Pass bulletin  as argument
$bull_date=$ARGV[0];
$bull_date=~s/f//;


# Dimensions -- lengths
%DimLen=(
   'ForecastLength'     , 10,
   'Region'   , 22,
   'MAXCHAR'  , 80,
   'Interval' , 2);

# Dimensions -- name
%DimName=(
   'ForecastLength'    ,'ForecastLength',
   'Region'  ,'Region',
   'MAXCHAR'  ,'MAXCHAR',
   'Interval','Interval' );

# Coordinate Variables ++  name
%VarName=(
   'ForecastLength'         ,'ForecastLength',
   'Region'       ,'Region',
   'RegionName'   ,'RegionName',
   'RegionLon'    ,'RegionLon',
   'RegionLat'    ,'RegionLat',
   'Interval'     ,'Interval',
   'DepthInterval','DepthInterval');


# Field (RMS and in time BIAS) Variables -- name
%FLDVarName=(
   'RMSINNOV1SST' ,'RMS_INNOV1_SST',
   'RMSINNOV1TEMP','RMS_INNOV1_TEMP',
   'RMSINNOV1SALN','RMS_INNOV1_SALN',
   'RMSINNOV1TIDE','RMS_INNOV1_TIDE',
   'RMSPERSSST'   ,'RMS_PERS_SST',
   'RMSPERSTEMP'  ,'RMS_PERS_TEMP',
   'RMSPERSSALN'  ,'RMS_PERS_SALN',
   'RMSPERSTIDE'  ,'RMS_PERS_TIDE',
   'RMSINNOVSST'  ,'RMS_INNOV_SST', 
   'RMSINNOVSLA'  ,'RMS_INNOV_SLA', 
   'RMSINNOVICEC' ,'RMS_INNOV_ICEC' ,
   'RMSRESIDSST'  ,'RMS_RESID_SST', 
   'RMSRESIDSLA'  ,'RMS_RESID_SLA', 
   'RMSRESIDICEC' ,'RMS_RESID_ICEC' );

# Variables -- associated file variable names (element of file name)
%FileVarName=(
   'RMSINNOV1SST' ,'SST',
   'RMSINNOV1TEMP','TEMP',
   'RMSINNOV1SALN','SALN',
   'RMSINNOV1TIDE','TIDE' ,
   'RMSPERSSST'   ,'SST',
   'RMSPERSTEMP'  ,'TEMP',
   'RMSPERSSALN'  ,'SALN',
   'RMSPERSTIDE'  ,'TIDE', 
   'RMSINNOVSST'  ,'SST', 
   'RMSINNOVSLA'  ,'SLA', 
   'RMSINNOVICEC' ,'icec',
   'RMSRESIDSST'  ,'SST', 
   'RMSRESIDSLA'  ,'SLA', 
   'RMSRESIDICEC' ,'icec' );


# RMS Variables -- associated RMS type names (element of file name)
%RmsTypeVarName=(
   'RMSINNOV1SST' ,'INNOV1',
   'RMSINNOV1TEMP','INNOV1',
   'RMSINNOV1SALN','INNOV1',
   'RMSINNOV1TIDE','INNOV1',
   'RMSPERSSST'   ,'PERSISTENCE',
   'RMSPERSTEMP'  ,'PERSISTENCE',
   'RMSPERSSALN'  ,'PERSISTENCE',
   'RMSPERSTIDE'  ,'PERSISTENCE',
   'RMSINNOVSST'  ,'INNOV', 
   'RMSINNOVSLA'  ,'INNOV', 
   'RMSINNOVICEC' ,'INNOV' ,
   'RMSRESIDSST'  ,'RESID', 
   'RMSRESIDSLA'  ,'RESID', 
   'RMSRESIDICEC' ,'RESID' );

%FLDVarUnit=(
   'RMSINNOV1SST' ,'degree_celius2',
   'RMSINNOV1TEMP','degree_celius2',
   'RMSINNOV1SALN','psu2',
   'RMSINNOV1TIDE','m2',
   'RMSPERSSST'   ,'degree_celius2',
   'RMSPERSTEMP'  ,'degree_celius2',
   'RMSPERSSALN'  ,'psu2',
   'RMSPERSTIDE'  ,'m2',
   'RMSINNOVSST'  ,'degree_celsius2', 
   'RMSINNOVSLA'  ,'m2', 
   'RMSINNOVICEC' ,'' ,
   'RMSRESIDSST'  ,'degree_celsius2', 
   'RMSRESIDSLA'  ,'m2', 
   'RMSRESIDICEC' ,'' );
 
%obsDim = (
   'RMSINNOV1SST'  , 2,
   'RMSINNOV1TIDE' , 2,
   'RMSINNOV1SALN' , 3,
   'RMSINNOV1TEMP' , 3,
   'RMSPERSSST'    , 2,
   'RMSPERSTIDE'   , 2,
   'RMSPERSSALN'   , 3,
   'RMSPERSTEMP'   , 3,
   'RMSINNOVSST'   , 1, 
   'RMSINNOVSLA'   , 1, 
   'RMSINNOVICEC'  , 1,
   'RMSRESIDSST'   , 1, 
   'RMSRESIDSLA'   , 1, 
   'RMSRESIDICEC'  , 1 );

%Average = (
   'RMSINNOV1SST'  , 1,
   'RMSINNOV1TIDE' , 1,
   'RMSINNOV1SALN' , 1,
   'RMSINNOV1TEMP' , 1,
   'RMSPERSSST'    , 1,
   'RMSPERSTIDE'   , 1,
   'RMSPERSSALN'   , 1,
   'RMSPERSTEMP'   , 1,
   'RMSINNOVSST'   , 0, 
   'RMSINNOVSLA'   , 0, 
   'RMSINNOVICEC'  , 0,
   'RMSRESIDSST'   , 0, 
   'RMSRESIDSLA'   , 0, 
   'RMSRESIDICEC'  , 0 );


# Create file name
$ncfile="topaz_V1_mersea_nat_rms_b${bull_date}.nc";

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
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "interpolation" , NetCDF::CHAR, 
   "Model values are interpolated to observation points");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "Regions" , NetCDF::CHAR, 
   "Observation points are collected in boxes defined by variables RegionLon, RegionLat");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "Depth_Interval" , NetCDF::CHAR, 
   "Profile data are averaged over depth intervals defined by variable DepthInterval");
$attid = NetCDF::attput($ncid, NetCDF::GLOBAL, "Time_Averaging" , NetCDF::CHAR, 
   "RMS Values are created by averaging over time interval [N-1 .. N+1] where N in [1,...,ForecastLength]");


# Create dimensions in NetCDF file
$DimId{'ForecastLength'}     = NetCDF::dimdef($ncid, $DimName{'ForecastLength'}    , $DimLen{'ForecastLength'} );
$DimId{'Region'}   = NetCDF::dimdef($ncid, $DimName{'Region'}  , $DimLen{'Region'} );
$DimId{'Interval'} = NetCDF::dimdef($ncid, $DimName{'Interval'}, $DimLen{'Interval'} );
$DimId{'MAXCHAR'}  = NetCDF::dimdef($ncid, $DimName{'MAXCHAR'} , $DimLen{'MAXCHAR'} );


# Define variables and attributes 
# 1) ForecastLength
$VarId{'ForecastLength'}= NetCDF::vardef($ncid, $VarName{'ForecastLength'}, NetCDF::FLOAT, $DimId{'ForecastLength'});
$attid = NetCDF::attput($ncid, $VarId{'ForecastLength'}, "long_name", NetCDF::CHAR, "Forecast Length in days ");
$attid = NetCDF::attput($ncid, $VarId{'ForecastLength'}, "axis", NetCDF::CHAR, "X");
$attid = NetCDF::attput($ncid, $VarId{'ForecastLength'}, "units", NetCDF::CHAR, "days");

# Define variables and attributes for
# 2) Region
@textdims=($DimId{'Region'},$DimId{'MAXCHAR'});
$VarId{'Region'}= NetCDF::vardef($ncid, $VarName{'Region'}, NetCDF::FLOAT, $DimId{'Region'});
$VarId{'RegionName'}= NetCDF::vardef($ncid, $VarName{'RegionName'}, NetCDF::CHAR, \@textdims);
$VarId{'RegionLon'}= NetCDF::vardef($ncid, $VarName{'RegionLon'}, NetCDF::CHAR, \@textdims);
$VarId{'RegionLat'}= NetCDF::vardef($ncid, $VarName{'RegionLat'}, NetCDF::CHAR, \@textdims);
$attid = NetCDF::attput($ncid, $VarId{'Region'}, "long_name", NetCDF::CHAR, "MERSEA Defined Region");
$attid = NetCDF::attput($ncid, $VarId{'Region'}, "axis", NetCDF::CHAR, "Y");
$attid = NetCDF::attput($ncid, $VarId{'RegionName'}, "long_name", NetCDF::CHAR, "Name of Region");
$attid = NetCDF::attput($ncid, $VarId{'RegionLon'}, "long_name", NetCDF::CHAR, "Region Boundary Longitude");
$attid = NetCDF::attput($ncid, $VarId{'RegionLat'}, "long_name", NetCDF::CHAR, "Region Boundary Latitude");


# Define variables and attributes for
# 3) Depth Interval
@textdims=($DimId{'Interval'},$DimId{'MAXCHAR'});
$VarId{'Interval'}= NetCDF::vardef($ncid, $VarName{'Interval'}, NetCDF::FLOAT, $DimId{'Interval'});
$VarId{'DepthInterval'}= NetCDF::vardef($ncid, $VarName{'DepthInterval'}, NetCDF::CHAR, \@textdims);
$attid = NetCDF::attput($ncid, $VarId{'Interval'}, "long_name", NetCDF::CHAR, "MERSEA Defined Interface");
$attid = NetCDF::attput($ncid, $VarId{'Interval'}, "axis", NetCDF::CHAR, "Z");
$attid = NetCDF::attput($ncid, $VarId{'DepthInterval'}, "long_name", NetCDF::CHAR, "Depth Interval Limits");



# RMS Variables -- name
@RMSDims3D=($DimId{'ForecastLength'},$DimId{'Region'},$DimId{'Interval'});
@RMSDims2D=($DimId{'ForecastLength'},$DimId{'Region'});
@RMSDims1D=($DimId{'Region'});
foreach $obstype (keys %FLDVarName) {
   if ($obsDim{$obstype}==3) {
      $FLDVarId{$obstype} = NetCDF::vardef($ncid, $FLDVarName{$obstype} , NetCDF::FLOAT, \@RMSDims3D);
   }elsif ($obsDim{$obstype} == 2) {
      $FLDVarId{$obstype} = NetCDF::vardef($ncid, $FLDVarName{$obstype} , NetCDF::FLOAT, \@RMSDims2D);
   }elsif ($obsDim{$obstype} == 1) {
      $FLDVarId{$obstype} = NetCDF::vardef($ncid, $FLDVarName{$obstype} , NetCDF::FLOAT, \@RMSDims1D);
   }
   $attid = NetCDF::attput($ncid, $FLDVarId{$obstype}, "_FillValue", NetCDF::FLOAT, $undef);
   $attid = NetCDF::attput($ncid, $FLDVarId{$obstype}, "units", NetCDF::CHAR, $FLDVarUnit{$obstype});
}



# End definition
$status = NetCDF::endef($ncid);



#Insert coordinate arrays
for ($region=1; $region<=$DimLen{'Region'} ; $region++) {
   $firstregion[$region]=1;
   $status = NetCDF::varput($ncid, $VarId{'Region'}, $region-1, 1, $region );
}
for ($intf=1; $intf<=$DimLen{'Interval'} ; $intf++) {
   $firstinterval[$intf]=1;
   $status = NetCDF::varput($ncid, $VarId{'Interval'}, $intf-1, 1, $intf );
}
for ($bin1=0; $bin1<$DimLen{'ForecastLength'} ; $bin1++) {
   $status = NetCDF::varput($ncid, $VarId{'ForecastLength'}, $bin1, 1, $bin1+1 );
}

# Gather data from RMS-files
foreach $obstype (keys %FLDVarId) {
   for ($region=1; $region<=$DimLen{'Region'} ; $region++) {



      $creg="00$region";
      $creg = substr($creg,length($creg)-3,3);
      $creg="r$creg";
      #print "$creg\n";

      if ($obsDim{$obstype} == 3) 
         { $intf_dim=$DimLen{'Interval'} }
      else 
         {$intf_dim=1}

      for ($intf=1; $intf<$intf_dim+1 ; $intf++) {
         $cint="00$intf";
         $cint = substr($cint,length($cint)-3,3);
         $cint="i$cint";
         #print "$cint\n";
         #print "$obstype $creg $cint $intf_dim\n";
         
         @infiles=();
         if ($obsDim{$obstype} == 3) {
            @infiles=glob("RMS.$RmsTypeVarName{$obstype}.b$bull_date" . 
                          "_f*.$FileVarName{$obstype}.$creg.$cint.dat"); 
         } else {
            @infiles=glob("RMS.$RmsTypeVarName{$obstype}.b$bull_date" . 
                          "_f*.$FileVarName{$obstype}.$creg.dat");
         }
         #print " Obs, 3dflag, int, reg, files:  $obstype $obsDim{$obstype} $cint $creg   @infiles\n";

         for ($bin1=0; $bin1 < $DimLen{'ForecastLength'} ; $bin1++) {
            $RMSvals[$bin1]=0;
            $RMSnrobs[$bin1]=0;
         }
         $maxbin=-1;

         # Loop through files
         for ($ifile=0; $ifile < $#infiles+1 ; $ifile++) {
            $filehandle = new IO::File;
            $filehandle->open("<$infiles[$ifile]") or die "Could not open file $infiles[$ifile]";

            # Skip header
            
            # Process Region name and index if first pass
            $line=$filehandle->getline();
            if ($firstregion[$region]==1) {
               $line=~ s/Region ://;
               @start=($region-1,0);
               $subline=substr($line,0,length($line));
               @count=(1,length($subline)-1);
               $status = NetCDF::varput($ncid, $VarId{'RegionName'}, \@start, \@count, $subline );
            }

            $filehandle->getline();
            $line=$filehandle->getline();
            # Process Longitude boundaries if first pass
            if ($firstregion[$region]==1) {
               $line=~ s/Region - Longitude ://;
               #print "$line \n";
               $subline=substr($line,0,length($line));
               @count=(1,length($subline)-1);
               $status = NetCDF::varput($ncid, $VarId{'RegionLon'}, \@start, \@count, $subline );
            }

            # Process Latitude boundaries if first pass
            $line=$filehandle->getline();
            if ($firstregion[$region]==1) {
               $line=~ s/Region - Latitude  ://;
               #print "$line \n";
               @count=(1,length($line)-1);
               $status = NetCDF::varput($ncid, $VarId{'RegionLat'}, \@start, \@count, $line );
               $firstregion[$region]=0;
            }


            # Process interval boundaries if first pass
            $line=$filehandle->getline();
            if ($firstinterval[$intf]==1 && $obsDim{$obstype}==3) {
               print "$line \n";
               $line=~ s/Depth interval ://;
               print "$line \n";
               @start=($intf-1,0);
               @count=(1,length($line)-1);
               $status = NetCDF::varput($ncid, $VarId{'DepthInterval'}, \@start, \@count, $line );
               $firstinterval[$intf]=0;
            }

            $filehandle->getline();

            # Actual data
            $line=$filehandle->getline();
            @flds=split " ","$line";
            $nrobs=$flds[6];
            $floattime=$flds[2];
            if (substr($obstype,0,3) eq "RMS") {
               $rmsval=$flds[5]; 
            }elsif (substr($obstype,0,4) eq "BIAS") {
               $rmsval=$flds[4]; 
            }else{
               print "Unknown observation type ...";
               exit
            }
            $filehandle->close;
            $bull_date=$flds[0];
            $forc_date=$flds[1];

            # Get bulletin and forecast day, month etc
            $forc_dim=substr($forc_date,6, 2);
            $bull_dim=substr($bull_date,6, 2);
            $forc_mnd=substr($forc_date,4, 2);
            $bull_mnd=substr($bull_date,4, 2);
            $forc_yr =substr($forc_date,0, 4);
            $bull_yr =substr($bull_date,0, 4);
            $delta=Delta_Days($bull_yr, $bull_mnd, $bull_dim, $forc_yr, $forc_mnd, $forc_dim) ;

            if ($verbose==1){
               print " $infiles[$ifile] $delta $FileVarName{$obstype} $region $intf $ifile $floattime $rmsval \n";
            }
            #print " $FLDVarId{$obstype} \n";
            #

            # Put into Array -- only forecast
            if ($delta>=0 && $obsDim{$obstype} > 1) {

               # This puts data into bins centered around [T0+1 T0+2 ...]
               # with range +-1. Or:  [T0+N-1, T0+N+1]
               #$delta=$delta-1;
               for ($bin1=0;$bin1<=2; $bin1++) {
                  $bin2=$bin1+$delta-1;
                  if ($bin2 < $DimLen{'ForecastLength'} && $bin2>=0 ) {

                     if (substr($obstype,0,3) eq "RMS") {
                        $RMSvals [$bin2]=$RMSvals [$bin2]+$rmsval**2*$nrobs;
                     }elsif (substr($obstype,0,4) eq "BIAS") {
                        $RMSvals [$bin2]=$RMSvals [$bin2]+$rmsval*$nrobs;
                     }
                     $RMSnrobs[$bin2]=$RMSnrobs[$bin2]+$nrobs;
                  }
               }


               if ($bin2 > $maxbin) {$maxbin=$bin2}
               #print "RMS nrobs: $RMSnrobs[$bin1]\n";
               if ($verbose==1) {
                  print " Delta bin1/2, dimeln, maxbin , file : $delta $bin1  $bin2 $DimLen{'ForecastLength'} $maxbin ";
                  print  " $infiles[$ifile]\n";
			   }
			} elsif ($delta==0 && $obsDim{$obstype}==1){
			   # Put straight into netcdf file...
			   @start=($region-1);
			   @count=(1);
			   if ($nrobs==0){ $rmsval=$undef;}
			   $status = NetCDF::varput($ncid, $FLDVarId{$obstype}, \@start, \@count, $rmsval );
			   print "$rmsval $nrobs $obstype\n";
			}
         }

         #print "RMS nrobs: @RMSnrobs\n";
         # Wrap up for this obs, region and interface
         # Put into NetCDF file
		 if ($obsDim{$obstype}>1){
			for ($bin1=0; $bin1 < $DimLen{'ForecastLength'} ; $bin1++) {
			   if ($obsDim{$obstype} == 3) {
				  @start=($bin1,$region-1,$intf-1);
				  @count=(1,1,1);
			   } else {
				  @start=($bin1,$region-1);
				  @count=(1,1);
			   }
			   if ($RMSnrobs[$bin1]==0){
				  $rmsval=$undef;
				  $status = NetCDF::varput($ncid, $FLDVarId{$obstype}, \@start, \@count, $rmsval );
			   } else{
				  if (substr($obstype,0,3) eq "RMS") {
					 $rmsval =sqrt($RMSvals[$bin1] / $RMSnrobs[$bin1]);
				  }elsif (substr($obstype,0,4) eq "BIAS") {
					 $rmsval =$RMSvals[$bin1] / $RMSnrobs[$bin1];
				  }
				  $status = NetCDF::varput($ncid, $FLDVarId{$obstype}, \@start, \@count, $rmsval );
			   }
			   #if ($verbose==1){ print "Dumping bins bin, nrobs, rms:  $bin1 $RMSnrobs[$bin1] $rmsval\n";}
			}
		 }


      }
   }
}

NetCDF::close($ncid) == 0 || die "Couldn't close netCDF file\n";
print "Fields dumped to file $ncfile\n";
