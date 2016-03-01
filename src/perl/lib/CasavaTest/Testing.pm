package CasavaTest::Testing;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw(&getTestXmlPath &processProjectTestXml &processCampaign 
      &processSuite &processTestSpecification &readPipelineTestXml 
      &readPipelineTestXml &getPipelineTestXmlpath &runTest 
      $ENV_NAME_GA_TEST_DATA $pipelineTestXmlName );
    @EXPORT_OK = qw();
}

# PROJECT: Pipeline
# MODULE:  Testing.pm
# AUTHOR:  L. Szajkowski
#
# Copyright (c)  2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Functions used to test and validate the pipeline
# Procedures in the library read the configuration and enviroment
# variables, they sepups the data and then delegates testing to
# runValidate, runGoat, runBustard or runGerald
#

=pod

=head1 NAME

CasavaTest::Testing.pm - Functions used to test and validate the pipeline

=head2 SYNOPSIS

use CasavaTest::Testing.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2007 Solexa, 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 DESCRIPTION

=head2 Overview

Messaging Functions 

=head2 Exports

Global variables:

=head2 Depends on

    warnings, strict, POSIX, XML::Simple, constant, Carp

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use XML::Simple;
use constant DEBUG => 0;    # set to 1 to get debug info
use Data::Dumper;
use Carp;

use File::Copy;
use CasavaTest::CasavaTesting;
#use CasavaTest::GeraldTesting;
#use CasavaTest::BustardTesting;
#use CasavaTest::FirecrestTesting;
#use CasavaTest::Validate;
use CasavaTest::TestingTools;

sub getTestXmlPath($;$);
sub processProjectTestXml($;$;$;$);
sub processCampaign($;$;$);
sub processSuite($;$;$;$);
sub processTestSpecification($;$;$;$);
sub runTest($;$;$;$;$;$);
sub readPipelineTestXml($);
sub readProjectTestXmlFile($);
sub getPipelineTestXmlpath($);
sub getGetSupportedTypes(\%;\%);

sub setUpTest(\%;$;$;$;$;$);
sub tearDownTest(\%;$;$;$;$;$);

our $ENV_NAME_GA_TEST_DATA = "GA_TEST_DATA";
our $ENV_NAME_GA_BIN_DIR   = "GA_BIN_DIR";
our $ENV_NAME_GA_DIR       = "GA_DIR";
our $pipelineTestXmlName   = "pipeline_tests.xml";
our $validateScript        = "validate.sh";

our $projectTestPrefix     = "ProjectTests-";
our $campaignPrefix        = "Campaign-";
our $suitePrefix           = "Suite-";
our $testPrefix            = "TestSpecification-";
our $elementContent        = "content";

# Project Test element id constants
our $testTestCampaign           = "TestCampaign";
our $testTestCampaignSuffix     = "-TestCampaigns";

# Test campaign element id constants
our $campaignCampaignName       = "campaignName";
our $campaignCampaignNameSuffix = "-$campaignCampaignName";
our $campaignDataFolder         = "runFolder";
our $campaignDataFolderSuffix   = "-$campaignDataFolder";
our $campaignTestSuite          = "testSuite";
our $campaignTestSuiteSuffix    = "-$campaignTestSuite";


# Test suite element id constants
our $suiteId                     = "id";
our $suiteIdSuffix               = "-$suiteId";
our $suiteCyclesDefault          = "cycles";
our $suiteCyclesDefaultSuffix    = "-$suiteCyclesDefault";
our $suiteCompatibility          = "compatibility";
our $suiteCompatibilitySuffix    = "-$suiteCompatibility";
our $suiteEnabled                = "enabled";
our $suiteEnabledSuffix          = "-$suiteEnabled";
our $suiteDataSetId              = "dataSet";
our $suiteDataSetIdSuffix        = "-suiteDataSetId";
our $suiteDataSetSuffix          = "-$suiteDataSetId";
our $suiteFasta                  = "fasta";
our $suiteFastaSuffix            = "-$suiteFasta";
our $suiteReadMode               = "readMode";
our $suiteReadModeSuffix         = "-$suiteReadMode";
our $suiteGerald                 = "GERALD";
our $suiteGeraldSuffix           = "-$suiteGerald";
our $suiteGeraldFolder           = "GERALDFolder";
our $suiteGeraldFolderSuffix     = "-$suiteGeraldFolder";
our $suiteType                   = "type";
our $suiteTypeSuffix             = "-$suiteType";
our $suiteTestSpecification      = "testSpecification";
our $suiteTestSpecificationSuffix= "-$suiteTestSpecification";

# Test specification element id contents
our $testspecId                  =  "id";
our $testspecIdSuffix            =  "-$testspecId";
our $testspecTag                 =  "tag";
our $testspecTagSuffix           =  "-$testspecTag";
our $testspecSet                 =  "set";
our $testspecSetSuffix           =  "-$testspecSet";
our $testspecCycles              =  "cycles";
our $testspecCyclesSuffix        =  "-$testspecCycles";
our $testspecLanes               =  "lanes";
our $testspecLanesSuffix         =  "-$testspecLanes";
our $testspecTiles               =  "tiles";
our $testspecTilesSuffix         =  "-$testspecTiles";
our $testspecEnabled             =  "enabled";
our $testspecEnabledSuffix       =  "-$testspecEnabled";
our $testspecScript              =  "script";
our $testspecScriptSuffix        =  "-$testspecScript";
our $testspecSetUp               =  "setUp";
our $testspecSetUpSuffix         =  "-$testspecSetUp";
our $testspecTearDown            =  "tearDown";
our $testspecTearDownSuffix      =  "-$testspecTearDown";
our $testspecSetUpPostCommand    =  "setUpPostCommand";
our $testspecSetUpPostCommandSuffix =  "-$testspecSetUpPostCommand";
our $testspecTearDownPostCommand =  "tearDownPostCommand";
our $testspecTearDownPostCommandSuffix =  "-$testspecTearDownPostCommand";
our $testspecTestCommand         =  "testCommand";
#our $testspecTestCommandSuffix   =  "-$testspecTestCommand";
our $testspecTestCommandSuffix   =  "-test_command";
our $testspecTestOptions         =  "commandOptions";
our $testspecTestOptionsSuffix   =  "-$testspecTestCommand";
our $testspecExpectedOutput      =  "expectedOutput";
our $testspecExpectedOutputSuffix=  "-$testspecExpectedOutput";
our $testspecExpectedError       =  "expectedError";
our $testspecExpectedErrorSuffix =  "-$testspecExpectedError";
our $testspecXpath               =  "xpath";
our $testspecXpathSuffix         =  "-$testspecXpath";
our $testspecTile                =  "tile";
our $testspecTileSuffix          =  "-testspecTile";
our $testspecLane                =  "lane";
our $testspecLaneSuffix          =  "-$testspecLane";
our $testspecAttributeCmd        =  "cmd";

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item info($msg)

The readPipelineTestXml reads pipeline_tests.xml

B<Parameters:>

    $path            - path to pipeline_tests.xml

B<Returns:> 

    HASH MAP Ref to struture with pipeline_tests.xml

=cut

sub readPipelineTestXml($) {
    croak("ERROR: readPipelineTestXml wrong number of parameters")
      unless ( @_ == 1 );
    my ($path) = @_;
    my $xs = new XML::Simple(
        searchpath => ".",
        forcearray => 1,
    );
    if ( !-e $path ) {
        die "ERROR: readPipelineTestXml file $path doesn't exist.\n";
    }
    my $ref = $xs->XMLin($path);
    return $ref;
}    # readPipelineTestXml

=pod

=item processPipelineTestXml($testsConfRef, $testId, $testSetId, $type)

The processPipelineTestXml puts default vaulues to test and test set 
definitions

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    HASH PAR Ref to test def, HASH PAR Ref to test set def  

=cut

sub processPipelineTestXml(\%;$;$;$) {
    croak("ERROR: processPipelineTestXml wrong number of parameters")
      unless ( @_ == 4 );
    my ( $testsConfRef, $testId, $testSetId, $type ) = @_;

    my $testRef = $testsConfRef->{test}->{$testId};
    if ( !defined $testRef ) {
        die "ERROR: $0 test [$testId] cannot be found in pipelineTestXml\n";
    }

    my $testSetRef = $testsConfRef->{set}->{$testSetId};
    if ( !defined $testSetRef ) {
        die "ERROR: $0 test set [$testSetId] cannot be found in pipelineTestXml\n";
    }

    my @geraldConfigs = @{ $testSetRef->{GERALD} };
    if ( scalar(@geraldConfigs) != 1 ) {
        die "ERROR: test framework supports one and only one gerald config for each data set\n";
    }

    my $typesRef = getGetSupportedTypes( %$testRef, %$testSetRef );
    my $isTypeSupported = 0;
    foreach my $supportedType (@$typesRef) {
        if ( $supportedType eq $type ) {
            $isTypeSupported = 1;
        }
    }

    if ( $isTypeSupported == 0 ) {
        die
"ERROR: $0 test type [$type] is not supported for $testId $testSetId\n";
    }

    my %testDefaults = %{ $testsConfRef->{conf}[0]->{testDefaults}[0] };

    #print Dumper( \%testDefaults );
    foreach my $key ( keys %testDefaults ) {
        if ( !defined $testRef->{$key} ) {

            #print "$key $testDefaults{$key}\n";
            $testRef->{$key} = $testDefaults{$key};
        }
    }

    my %testSetDefaults = %{ $testsConfRef->{conf}[0]->{testSetDefaults}[0] };
    foreach my $key ( keys %testSetDefaults ) {
        if ( !defined $testSetRef->{$key} ) {
            $testSetRef->{$key} = $testSetDefaults{$key};
        }
    }

    #print Dumper($testSetRef);
    return $testRef, $testSetRef;
}    # processPipelineTestXml

=pod

=item getPipelineTestXmlpath($msg)

The getPipelineTestXmlpath generates path to pipeline_tests.xml file
GA_TEST_DATA can be set either from command line or in enviroment.

B<Parameters:>

    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)

B<Returns:> 

    path to pipeline_tests.xml file

=cut

sub getPipelineTestXmlpath($) {
    croak("ERROR: getPipelineTestXmlpath wrong number of parameters")
      unless ( @_ == 1 );
    my ($GA_TEST_DATA_PATH) = @_;
    if ( $GA_TEST_DATA_PATH eq "" && defined $ENV{$ENV_NAME_GA_TEST_DATA} ) {
        $GA_TEST_DATA_PATH = $ENV{$ENV_NAME_GA_TEST_DATA};
    }

#    print "The test data path is %GA_TEST_DATA_PATH\n";

    if ( !-d $GA_TEST_DATA_PATH ) {
        die
            "ERROR: $0 $ENV_NAME_GA_TEST_DATA=$GA_TEST_DATA_PATH is not accessible from getPipelineTestXmlpath()\n";
    }

    my $pipelineTestXmlPath =
      File::Spec->catfile( $GA_TEST_DATA_PATH, $pipelineTestXmlName );


    if ( !-e $pipelineTestXmlPath ) {
        die "ERROR: $0 the pipelineTest XML file $pipelineTestXmlPath is not accessible\n";
    }

    return ( $pipelineTestXmlPath, $GA_TEST_DATA_PATH );
}    # getPipelineTestXmlpath

=pod

=item
=cut

sub getTestXmlPath($;$)
{
    croak("ERROR: getTestXmlpath - wrong number of parameters")
      unless (@_ == 2);

    my ($GA_TEST_DATA_PATH, $testXmlFileName) = @_;
#    print "The project XML file name is $testXmlFileName\n";
    if ( $GA_TEST_DATA_PATH eq "" && defined $ENV{$ENV_NAME_GA_TEST_DATA} ) 
    {
        $GA_TEST_DATA_PATH = $ENV{$ENV_NAME_GA_TEST_DATA};
    }

    if ( !-d $GA_TEST_DATA_PATH ) {
        die
            "ERROR: $0 $ENV_NAME_GA_TEST_DATA=$GA_TEST_DATA_PATH is not accessible from getTestXmlPath()\n";
    }

    my $testXmlPath = File::Spec->catfile($GA_TEST_DATA_PATH, "testXml/".$testXmlFileName);

    if(! -e $testXmlPath)
    {
        die "ERROR: $0 the project test XML file $ENV_NAME_GA_TEST_DATA=$GA_TEST_DATA_PATH with full path $testXmlPath is not accessible\n";
    }

    return ($testXmlPath, $GA_TEST_DATA_PATH);
}

=pod

=item readProjectTestXmlFile()

Reads a project test configuration file specified by $path

B<Parameters:>

    $path            - path to a project test xml file

B<Returns:>

    HASH MAP Ref to structure containing a project tests xml tree

=cut

sub readProjectTestXmlFile($)
{
    croak("ERROR: readProjectTestXmlFile - wrong number of parameters")
      unless (@_ == 1);
    my ($path) = @_;
    my $xs = new XML::Simple(
        searchpath => ".",
        forcearray => 1,
    );
    if (!-e $path) {
        die "ERROR: readProjectTestXmlFile - file $path doesn't exist.\n";
    }
    my $ref = $xs->XMLin($path);
    return $ref;
}

=pod

=item processConfigurationXmlFile()

Extracts default configuration parameters from the XML file specified by $path
and places them in an associative array for later processing.

B<Parameters:>

    %testRef         - Hash map reference to structure containing test parameters
    $path            - path to a project test xml file

B<Returns:>

    HASH MAP Ref to structure containing a project tests xml tree

=cut

sub processConfigurationXmlFile($;$)
{
    croak("ERROR: processConfigurationXmlFile - wrong number of parameters")
      unless (@_ == 2);
    my ($testRef, $path) = @_;

    my $defConRef = readPipelineTestXml($path);
    if(!defined $defConRef)
    {
        die "ERROR: processConfigurationXmlFile - file $path could not be read.\n";
    }
    
    my $defaultConfigValues;

    foreach my $elementKey(keys %{$defConRef})
    {
        my %testDefaults = %{$defConRef->{$elementKey}[0]};
        foreach my $attributeKey(keys %testDefaults)
        {
            my $value = $testDefaults{$attributeKey};
            if(!defined $defaultConfigValues->{$elementKey}{$attributeKey})
            {
                $defaultConfigValues->{$elementKey}{$attributeKey} = $value;
            }
        }
    }
    return $defaultConfigValues;
}

=pod
=item
   <WEB_DIR_ROOT address="http://snark/" />
   <errorReport tearDown="Y" noveData="0" />
   <testDefaults cycles="C1.1,C2.1,C3.1,C4.1,C5.1,C6.1,C7.1,C8.1" 
                 lanes="1,2,3" 
                 tiles="1,2,3" 
                 runMake="N"
                 setUp="Y"
                 tearDown="Y"
                 compatibility="1.1rc1p4"
                 enabled="Y"
                 addCmdOption="" />
   <testSetDefaults enabled="Y"> 
=pod

=item processProjectTestXml($testsConfRef, $testId, $testSetId, $type)

The processProjectTestXml puts default vaulues to test and test set 
definitions

B<Parameters:>

    $testRef           - HASH MAP Ref to struture containing <Project>Test.xml
    $path              - path to the <project>Test.xml file
    $testSetId         - Name of the test set to be processed
    $testParameterData - reference to a HASH structure used to store parameter details related to the test run 

B<Returns:> 

    ARRAY Ref to test campaigns specified in the project test file

=cut

sub processProjectTestXml($;$;$;$) {
    croak("ERROR: processProjectTestXml wrong number of parameters")
      unless ( @_ == 4 );
    my ($testRef, $path, $testSetId, $testParameterData) = @_;
    #print "In processProjectTestXml\n";

    my $campaignsRef = readPipelineTestXml($path);
    if(!defined $campaignsRef)
    {
        die "ERROR: processProjectTestXml - file $path could not be read.\n";
    }
    
    my $campaignList;

    foreach my $campaign(@{$campaignsRef->{$testTestCampaign}})
    {
#        print "Campaign name - $campaign\n";
        push(@$campaignList, $campaign);
    }
    $testParameterData->{$projectTestPrefix.$testSetId.$testTestCampaignSuffix} = $campaignList;

    return $campaignList, $testParameterData;
}

=pod

=item processTestCampaignXml($testsConfRef, $testId, $testSetId, $type)

The processTestCampaignXml puts default vaulues to test and test set 
definitions

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    HASH PAR Ref to test def, HASH PAR Ref to test set def  

=cut

sub processTestCampaignXml($;$;$) {
    #print "In processTestCampaignXml()\n";
    croak("ERROR: processTestCampaignXml wrong number of parameters")
      unless ( @_ == 3 );
    my ($campaign, $path, $testSpecificationData ) = @_;

    my $campaignsRef = readPipelineTestXml($path);
    if(!defined $campaignsRef)
    {
        die "ERROR: processTestCampaignXml - file $path could not be read.\n";
    }

    my $campaignName = $campaignsRef->{$campaignCampaignName};
    if(defined $campaignName)
    {
        $testSpecificationData->{$campaignPrefix.$campaign.$campaignCampaignNameSuffix} = $campaignName;
    }

    foreach my $dataFolderName (@{$campaignsRef->{$campaignDataFolder}})
    {
        my $folderKey = $campaignPrefix.$campaign.$campaignDataFolderSuffix;
        $testSpecificationData->{$folderKey} = $dataFolderName;
#        print "Element $folderKey contains the data folder name - $dataFolderName\n";
    }
    my $testSuiteList;
    foreach my $suite(@{$campaignsRef->{$campaignTestSuite}})
    {
#        print "Test suite name - $suite\n";
        push(@{$testSuiteList}, $suite);
    }
    $testSpecificationData->{$campaignPrefix.$campaign.$campaignTestSuiteSuffix} = $testSuiteList;

    return $testSuiteList;
}

=pod

=item processTestSuiteXml($testsConfRef, $testId, $testSetId, $type)

The processTestSuiteXml puts default vaulues to test and test set 
definitions

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    HASH PAR Ref to test def, HASH PAR Ref to test set def  

=cut

sub processTestSuiteXml($;$;$) {
    croak("ERROR: processTestSuiteXml wrong number of parameters")
      unless (@_ == 3);
    my ($suite, $testParameterData, $path) = @_;

    my $suiteRef = readPipelineTestXml($path);
    if(!defined $suiteRef)
    {
        die "ERROR: processTestSuiteXml - file $path could not be read.\n";
    }

    my $testSuiteId = $suiteRef->{$suiteId};
    if(defined $testSuiteId)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteIdSuffix} = $testSuiteId;
    }

    my $defaultCycles = $suiteRef->{$suiteCyclesDefault};
    if(defined $defaultCycles)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteCyclesDefaultSuffix} = $defaultCycles;
    }

    my $dataSetElement = $suiteRef->{$suiteDataSetId};
    if(defined $dataSetElement)
    {
        foreach my $compat(keys %$dataSetElement)
        {
            foreach my $element(keys %{$dataSetElement->{$compat}})
            {
                my $value = $dataSetElement->{$compat}->{$element};
                if(length($value) > 0)
                {
                    $testParameterData->{$suitePrefix.$suite."-".$element} = $value;
                }
                else
                {
                    $testParameterData->{$suitePrefix.$suite."-".$element} = undef;
                }
            }
        }
    }

    my $dataSetCompatibility = $suiteRef->{$suiteDataSetId.$suiteCompatibility};
    if(defined $dataSetCompatibility)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteCompatibilitySuffix} = $dataSetCompatibility;
    }

    my $isEnabled = $suiteRef->{$suiteEnabled};
    if(defined $isEnabled)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteEnabledSuffix} = $isEnabled;
    }
    else
    {
        print "Test suite has not been specified as enabled or not, has been set to default of NOT enabled\n";
        $testParameterData->{$suitePrefix.$suite.$suiteEnabledSuffix} = "N";
    }

    my $datasetId = $suiteRef->{$suiteDataSetId};
    my $datasets;
    if(defined $datasetId)
    {
        foreach my $set(keys %{$datasetId})
        {
            my $dataset = $datasetId->{$set}->{$elementContent};
            $datasets->{$set} = $dataset;
        }
        $testParameterData->{$suitePrefix.$suite.$suiteDataSetSuffix} = $datasets;
    }

    my $fastaFile = $suiteRef->{$suiteFasta};
    my $fastaFileList;
    if(defined $fastaFile)
    {
        foreach my $file(@{$fastaFile})
        {
            push(@{$fastaFileList}, $file);
        }
        $testParameterData->{$suitePrefix.$suite.$suiteFastaSuffix} = $fastaFileList;
    }

    my $readMode = $suiteRef->{$suiteReadMode};
    if(defined $readMode)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteReadModeSuffix} = $readMode;
    }

    my $gerald = $suiteRef->{$suiteGerald};
    my $geraldFileList;
    if(defined $gerald)
    {
        foreach my $geraldFile(@{$gerald})
        {
            push(@{$geraldFileList}, $geraldFile);
#            print "\t\tGERALD config file - $geraldFile\n";
        }
        $testParameterData->{$suitePrefix.$suite.$suiteGeraldSuffix} = $geraldFileList;
    }

    my $geraldFolder = $suiteRef->{$suiteGeraldFolder};
    if(defined $geraldFolder)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteGeraldFolderSuffix} = $geraldFolder;
    }

    my $type = $suiteRef->{$suiteType};
    if(defined $type && length($type->[0]) > 0)
    {
        $testParameterData->{$suitePrefix.$suite.$suiteTypeSuffix} = $type;
    }
    else
    {
        push(@$type, "Not Specified");
        $testParameterData->{$suitePrefix.$suite.$suiteTypeSuffix} = $type;
    }

    my $testSpecificationList;
    foreach my $testSpec(@{$suiteRef->{$suiteTestSpecification}})
    {
        push(@{$testSpecificationList}, $testSpec);
    }
    $testParameterData->{$suitePrefix.$suite.$suiteTestSpecificationSuffix} = $testSpecificationList;

    return $testSpecificationList;
}

=pod

=item processTestSpecificationXml($testsConfRef, $testId, $testSetId, $type)

The processTestSpecificationXml puts default vaulues to test and test set 
definitions

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    HASH PAR Ref to test def, HASH PAR Ref to test set def  

=cut

sub processTestSpecificationXml($;$;$) {
    croak("ERROR: processTestSpecificationXml wrong number of parameters")
      unless ( @_ == 3 );
    my ( $name, $testParameterData, $path ) = @_;

    my $testParamRef = readPipelineTestXml($path);
    if(!defined $testParamRef)
    {
        die "ERROR: processTestSpecificationXml - file $path could not be read.\n";
    }

    my $testParameterList;
    my $testParameters;
    my $specificationName = "$testPrefix$name";

    my $testId = $testParamRef->{$testspecId};
    if(defined $testId)
    {
        $testParameterData->{$specificationName.$testspecIdSuffix} = $testId;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecIdSuffix} = "NoTestId";
    }

    my $testTag = $testParamRef->{$testspecTag};
    if(defined $testTag)
    {
        $testParameterData->{$specificationName.$testspecTagSuffix} = $testTag;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecTagSuffix} = " ";
    }

    my $testSet = $testParamRef->{$testspecSet};
    if(defined $testSet)
    {
        $testParameterData->{$specificationName.$testspecSetSuffix} = $testSet;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecSetSuffix} = "N";
    }

    my $testCycles = $testParamRef->{$testspecCycles};
    if(defined $testCycles)
    {
        $testParameterData->{$specificationName.$testspecCyclesSuffix} = $testCycles;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecCyclesSuffix} = "N";
    }

    my $testLanes = $testParamRef->{$testspecLanes};
    if(defined $testLanes)
    {
        $testParameterData->{$specificationName.$testspecLanesSuffix} = $testLanes;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecLanesSuffix} = "N";
    }

    my $testTiles = $testParamRef->{$testspecTiles};
    if(defined $testTiles)
    {
        $testParameterData->{$specificationName.$testspecTilesSuffix} = $testTiles;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecTilesSuffix} = "N";
    }

    my $testSetUp = $testParamRef->{$testspecSetUp};
    if(defined $testSetUp)
    {
        $testParameterData->{$specificationName.$testspecSetUpSuffix} = $testSetUp;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecSetUpSuffix} = "N";
    }

    my $testSetUpPostCommand = $testParamRef->{$testspecSetUpPostCommand}[0]->{$testspecAttributeCmd};
    if(defined $testSetUpPostCommand)
    {
        $testParameterData->{$specificationName.$testspecSetUpPostCommandSuffix} = $testSetUpPostCommand;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecSetUpPostCommandSuffix} = "";
    }

    my $testTearDown = $testParamRef->{$testspecTearDown};
    if(defined $testTearDown)
    {
        $testParameterData->{$specificationName.$testspecTearDownSuffix} = $testTearDown;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecTearDownSuffix} = "N";
    }

    my $testTearDownPostCommand = $testParamRef->{$testspecTearDownPostCommand}[0]->{$testspecAttributeCmd};
    if(defined $testTearDownPostCommand)
    {
        $testParameterData->{$specificationName.$testspecTearDownPostCommandSuffix} = $testTearDownPostCommand;
    }
    else
    {
        $testParameterData->{$specificationName.$testspecTearDownPostCommandSuffix} = "";
    }

    my $testEnabled = $testParamRef->{$testspecEnabled};
    if(defined $testEnabled)
    {
        $testParameterData->{$specificationName.$testspecEnabledSuffix} = $testEnabled;
    }
    else
    {
#        print "Test specification has not been specified as enabled or not, has been set to default of NOT enabled\n";
        $testParameterData->{$specificationName.$testspecEnabledSuffix} = "N";
    }

    my $testScript = $testParamRef->{$testspecScript};
    if(defined $testScript)
    {
        $testParameterData->{$specificationName.$testspecScriptSuffix} = $testScript;
    }
    else
    {
#        print "\nCannot perform test - A test script is not specified\n"; 
    }

#    foreach my $testSpec(@{$testParamRef->{$testspecTestCommand}})
#    {
#        push(@{$testParameterList}, $testSpec);
#    }
#    $testParameterData->{$specificationName.$testspecTestCommandSuffix} = $testParamRef->{$testspecTestCommand};
#    @$testParameters{$testspecTestCommand} = $testParameterList;

    my $testCommand = $testParamRef->{$testspecTestCommand}[0]->{$testspecAttributeCmd};
    if(defined $testCommand)
    {
        $testParameterData->{$specificationName.$testspecTestCommandSuffix} = $testCommand;
#        print "The defined test command is - $testCommand\n";
    }
    else
    {
        $testParameterData->{$specificationName.$testspecTestCommandSuffix} = "";
#        print "The defined test command is - undefined\n";
    }

    my $testOptionList;
    foreach my $testSpec(@{$testParamRef->{$testspecTestOptions}})
    {
        push(@{$testOptionList}, $testSpec);
    }
    $testParameterData->{$specificationName.$testspecTestOptionsSuffix} = $testOptionList;
    @$testParameters{$testspecTestOptions} = $testOptionList;

    my $testOutputList;
    foreach my $testSpec(@{$testParamRef->{$testspecExpectedOutput}})
    {
        push(@{$testOutputList}, $testSpec);
    }
    $testParameterData->{$specificationName.$testspecExpectedOutputSuffix} = $testOutputList;
    @$testParameters{$testspecExpectedOutput} = $testOutputList;

    my $testErrorList;
    foreach my $testSpec(@{$testParamRef->{$testspecExpectedError}})
    {
        push(@{$testErrorList}, $testSpec);
    }
    $testParameterData->{$specificationName.$testspecExpectedErrorSuffix} = $testErrorList;
    @$testParameters{$testspecExpectedError} = $testErrorList;

    return $testParameterData;
}

=pod

=item getGetSupportedTypes($testRef, $testSetRef)

The getGetSupportedTypes merges types from test and test set and returns 
supported types.

B<Parameters:>

    $testRef     - HASH MAP Ref to test
    $testSetRef  - HASH MAP Ref to test set

B<Returns:> 

    ARRAY Ref to supported types

=cut

sub getGetSupportedTypes(\%;\%) {
    croak("ERROR: getGetSupportedTypes wrong number of parameters")
      unless ( @_ == 2 );
    my ( $testRef, $testSetRef ) = @_;
    my @testSetTypes = @{ $testSetRef->{type} };
    my @types        = ();
    foreach my $type ( @{ $testRef->{type} } ) {
        foreach my $testSetType (@testSetTypes) {
            if ( $type eq $testSetType ) {
                push @types, $testSetType;
            }
        }
    }
    return \@types;
}    # getGetSupportedTypes

=pod

=item setUpTest($GA_TEST_DATA_PATH, $testId, $testSetId, $type, $testRef, 
    $testSetRef)

The procedure set ups data required to run a test (copy the data) 
used by validate.sh script

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)
    $ENV_NAME_GA_DIR_PATH  - Path to pipeline installation
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    Status (0 | -1)

=cut

sub setUpTest(\%;$;$;$;$;$) {
    croak( "ERROR: setUpTest wrong number of parameters " . scalar(@_) . "\n" )
      unless ( @_ == 6 );
    my ( $testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH, $testId,
        $testSetId, $type )
      = @_;

    print "Content of testsConfRef -\n";
    print Dumper($testsConfRef);

    my $testRef    = $testsConfRef->{test}->{$testId};
    my $testSetRef = $testsConfRef->{set}->{$testSetId};

    print "Content of testSetRef -\n";
    print Dumper($testSetRef);

    print "Content of testRef -\n";
    print Dumper($testRef);

    my @runFolders = @{ $testSetRef->{runFolder} };
    ## test framework supports one and only one gerald config for each data set
    my $geraldConfig      = @{ $testSetRef->{GERALD} }[0];
    my $isIpar            = 0;
    my $testSetFolderPath =
      File::Spec->catdir( $ENV_NAME_GA_DIR_PATH, $tempSpace, $testSetId );
    my ( $lanesStr, $cyclesStr, $tilesStr );

    $lanesStr  = $testRef->{lanes};
    $tilesStr  = $testRef->{tiles};
    $cyclesStr = $testRef->{cycles};

    if ( !-d $testSetFolderPath ) {
        mkdir $testSetFolderPath;
    }
    if ( $type eq $IPAR_DIR_NAME ) {
        $isIpar = 1;
    }
    my %dirs = ();

    ## Copy images and runFolder basic structure
    my $runfoldersStr = "";
    my $gerald        = File::Spec->catdir( $runFolders[0], $geraldConfig );

    if ( $type eq $IPAR_DIR_NAME ) {
        if ( scalar(@runFolders) != 1 ) {
            die "ERROR: $0 $type run can have only one run folder\n";
        }
        my $runDir =
          File::Spec->catdir( $GA_TEST_DATA_PATH, $testSetId, $runFolders[0] );
        copyImageAnaFiles( $runDir, $testSetFolderPath, 1, %dirs );
        if ( $testRef->{script} eq "goat" || $testRef->{script} eq "bustard") {
            copyIpar( $runDir, $lanesStr, $cyclesStr, $tilesStr, %dirs );
        }
    }
    elsif ( $type eq "GA" ) {
        foreach my $runFolder (@runFolders) {
            my $runDir =
              File::Spec->catdir( $GA_TEST_DATA_PATH, $testSetId, $runFolder );
            copyImageAnaFiles( $runDir, $testSetFolderPath, 0, %dirs );
            if ( $testRef->{script} eq "goat" ) {
                copyImages( $runDir, $lanesStr, $cyclesStr, $tilesStr, %dirs );
            }
        }
    }
    else {
        die "ERROR: $0 Unsuported type [$type] in $testId test\n";
    }

    ## Set upt PairedEndInfo.xml only for SFPE runs
    if (   $cyclesStr ne "auto"
        && scalar(@runFolders) == 1
        && $testSetRef->{readMode}[0] eq "PE" )
    {
        my @cyclesArray = split ",", $cyclesStr;
        my $read1Length = int( scalar(@cyclesArray) / 2 );
        my $content     =
          "<?xml version=\"1.0\"?>\n<FirstRead Length=\"$read1Length\" />\n";
        my $pairedEndInfoPath =
          File::Spec->catfile( $testSetFolderPath, $runFolders[0],
            "PairedEndInfo.xml" );
        write2file( $pairedEndInfoPath, $content );
    }

    if ( $testRef->{script} ne "goat" ) {
        copyFirecrestAll( %$testsConfRef, $GA_TEST_DATA_PATH,
            $ENV_NAME_GA_DIR_PATH, $testId, $testSetId, $type );
    }

    my $geraldFolderRef = getGeraldFolder( $type, %$testRef, %$testSetRef );
    my $dirsRef         = parseGeraldPath( $geraldFolderRef->{path} );

    my $topFolderPath;

    my $script = $testRef->{script};
    if ( $script eq "validate" || $script eq "goat" ) {

    }
    elsif ( $script eq "bustard" ) {
        $topFolderPath =
          File::Spec->catdir( $testSetFolderPath, $dirsRef->{FirecrestFull} );
    }
    elsif ( $script eq "GERALD" ) {
        $topFolderPath =
          File::Spec->catdir( $testSetFolderPath, $dirsRef->{BustardFull} );
    }
    else {
        die "ERROR: unsupported script in [$testId] test\n";
    }

    if ( defined $testRef->{setUpPostCmd} ) {
        my @setUpPostCmd = @{ $testRef->{setUpPostCmd} };
        my $curCir       = Cwd->getcwd();
        if ( $script eq "validate" || $script eq "goat" ) {
            foreach my $runFolder (@runFolders) {
                chdir File::Spec->catdir( $testSetFolderPath, $runFolder );
                executeCmds(@setUpPostCmd);
            }
        }
        else {
            chdir $topFolderPath;
            executeCmds(@setUpPostCmd);
        }
        chdir File::Spec->catdir($curCir);
    }
    return 0;
}    # sub setUpTest

=pod

=item tearDownTest($GA_TEST_DATA_PATH, $testId, $testSetId, $type, $testRef, 
    $testSetRef)

The procedure cleans data after running a test

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)
    $ENV_NAME_GA_DIR_PATH  - Path to pipeline installati" . scalar(@_) . "\n"on
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    Status (0 | -1)

=cut

sub tearDownTest(\%;$;$;$;$;$) {
    croak(
        "ERROR: tearDownTest wrong number of parameters " . scalar(@_) . "\n" )
      unless ( @_ == 6 );
    my ( $testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH, $testId,
        $testSetId, $type )
      = @_;

    my $testRef           = $testsConfRef->{test}->{$testId};
    my $testSetRef        = $testsConfRef->{set}->{$testSetId};
    my @runFolders        = @{ $testSetRef->{runFolder} };
    my @fastas            = @{ $testSetRef->{fasta} };
    my $cycles            = $testRef->{cycles};
    my @geraldConfigs     = @{ $testSetRef->{GERALD} };      # <<<< Do we need anything for CASAVA configs???
    my $content           = "";
    my $testSetFolderPath =
      File::Spec->catdir( $ENV_NAME_GA_DIR_PATH, $tempSpace, $testSetId );

    if ( !-d $testSetFolderPath ) {
        ## There is no need to remove the data since each execution creates
        # a unique directory
        #executeCmd("rm -fr $testSetFolderPath");
    }

    my $geraldFolderRef = getGeraldFolder( $type, %$testRef, %$testSetRef );
    my $dirsRef = parseGeraldPath( $geraldFolderRef->{path} );
    my $topFolderPath;

    my $script = $testRef->{script};
    if ( $script eq "validate" ) {
    }
    elsif ( $script eq "goat" ) {
    }
    elsif ( $script eq "bustard" ) {
        $topFolderPath =
          File::Spec->catdir( $testSetFolderPath, $dirsRef->{FirecrestFull} );
    }
    elsif ( $script eq "GERALD" ) {
        $topFolderPath =
          File::Spec->catdir( $testSetFolderPath, $dirsRef->{BustardFull} );
    }
    elsif ( $script eq "Casava" ) 
    {
        $topFolderPath =
          File::Spec->catdir( $testSetFolderPath, $dirsRef->{CasavaFull} );
    }
    else {
        die "ERROR: unsupported script in [$testId] test\n";
    }

    if ( defined $testRef->{tearDownPostCmd} ) {    
        my $curCir          = Cwd->getcwd();
        my @tearDownPostCmd = @{ $testRef->{tearDownPostCmd} };
        if ( $script eq "validate" || $script eq "goat" ) {
            foreach my $runFolder (@runFolders) {
                chdir File::Spec->catdir( $testSetFolderPath, $runFolder );
                executeCmds(@tearDownPostCmd);
            }
        }
        else {
            chdir $topFolderPath;
            executeCmds(@tearDownPostCmd);
        }
    }
    return $content;
}    # sub tearDownTest

=pod

=item

=cut

sub processCampaign($;$;$)
{
    croak("ERROR: processCampaign wrong number of parameters")
      unless (@_ == 3);

    my ($testSpecificationData, $campaign, $GA_TEST_DATA_PATH) = @_;
#    print "Campaign to process - $campaign\n";

    my $thisCampaign;
    ($thisCampaign, $GA_TEST_DATA_PATH) = getTestXmlPath($GA_TEST_DATA_PATH, $campaign."_campaign.xml");
#    print "\nCampaign XML file is - $thisCampaign\n";
    
    my $testSuites = processTestCampaignXml($campaign, $thisCampaign, $testSpecificationData);
    my $folderKey = $campaignPrefix.$campaign.$campaignDataFolderSuffix;
    my $dataFolder;
    if(defined $testSpecificationData->{$folderKey})
    {
        $dataFolder = $testSpecificationData->{$folderKey};
#        print "\tData folder is - $dataFolder\n";
    }
    foreach my $key(keys %$testSpecificationData)
    {
#        print "Campaign parameter - $key\n";
        if($key =~ "testSuite")
        {
            my $suiteList = $testSpecificationData->{$key};
        }
    }
    return $testSuites;
}

=pod

=item
=cut

sub processSuite($;$;$;$)
{
    #print "In processSuite()\n";
    croak("ERROR: processSuite wrong number of parameters")
      unless (@_ == 4);

    my ($testSpecificationData, $suite, $dataFolder, $GA_TEST_DATA_PATH) = @_;

    my $thisSuite;
    ($thisSuite, $GA_TEST_DATA_PATH) = getTestXmlPath($GA_TEST_DATA_PATH, $suite.".xml");
    my $testSpecificationList = processTestSuiteXml($suite, $testSpecificationData, $thisSuite);

    if(!defined $testSpecificationList)
    {
        die "Test specification list was not available for suite $thisSuite\n";
    }

    return $testSpecificationList;
}

=pod

=item
=cut

sub processTestSpecification($;$;$;$)
{
    croak("ERROR: processTestSpecification wrong number of parameters")
      unless (@_ == 4);

    my ($testParameterData, $test, $dataFolder, $GA_TEST_DATA_PATH) = @_;

    my $thisTestSpecification;
    ($thisTestSpecification, $GA_TEST_DATA_PATH) = getTestXmlPath($GA_TEST_DATA_PATH, $test.".xml");
    my $specificationParameters = processTestSpecificationXml($test, $testParameterData, $thisTestSpecification);

    if(!defined $specificationParameters)
    {
        die "Test specification parameters are not available for test $thisTestSpecification\n";
    }

    return $specificationParameters; #$testParameterData;    
}

=pod

=item runTest($GA_TEST_DATA_PATH, $testId, $testSetId, $type)

The runTest runs test from pipeline_tests.xml

B<Parameters:>

    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    status (0, -1)

=cut

sub runTest($;$;$;$;$;$) {
    croak("ERROR: runTest wrong number of parameters")
      unless ( ( @_ == 3 ) || ( @_ == 4 ) || ( @_ == 6 ) );
    my ( $GA_TEST_DATA_PATH, $testId, $testSetId, $type, $testRef, $testSetRef ) = @_;
    my $pipelineTestXmlPath;
    my $testConfigXmlPath;
    my $testParameterData;
    my $status = -1;
    my $msg    = "";

    my $ENV_NAME_GA_DIR_PATH = "";
    if ( defined $ENV{$ENV_NAME_GA_DIR} ) {
        $ENV_NAME_GA_DIR_PATH = $ENV{$ENV_NAME_GA_DIR};
    }
    else {
        die "ERROR: Enviroment variable $ENV_NAME_GA_DIR not defined\n";
    }
        
    if($type eq "CASAVA")
    {
        print "Arguments supplied are\n\tGA_TEST_DATA_PATH = $GA_TEST_DATA_PATH\n";
        print "\tTest ID = $testId\n\tTest set ID = $testSetId\n\tTest type = $type\n";
        print Dumper($testRef);
        print Dumper($testSetRef);
        my $testReference;

        my $configFile;
        if(defined $testId."Configuration.xml")
        {
            $configFile = $testId."Configuration.xml"
        }
        else
        {
            $configFile = "defaultConfiguration.xml";
        }

        ($testConfigXmlPath, $GA_TEST_DATA_PATH) = getTestXmlPath($GA_TEST_DATA_PATH, $configFile);

        my $testsConfRef = processConfigurationXmlFile($testReference, $testConfigXmlPath);

        if ( $testRef->{enabled} eq "N" || $testSetRef->{enabled} eq "N" ) {
            print "test $testId for $type $testSetId disabled\n";
            return 0;
        }

        $testsConfRef->{test}->{$testId} = $testRef;
        $testsConfRef->{set}->{$testSetId} = $testSetRef;
        
        if ( $testRef->{setUp} eq "Y" ) {
            setUpTest( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                       $testId, $testSetId, $type );
        }
        
        ($status, $msg) = runCasava(%$testsConfRef, 
                                    $GA_TEST_DATA_PATH, 
                                    $ENV_NAME_GA_DIR_PATH,
                                    $testId, 
                                    $testSetId, 
                                    $type, 
                                    %$testRef, 
                                    %$testSetRef );

        if ( $testRef->{tearDown} eq "Y" ) {
            tearDownTest( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                          $testId, $testSetId, $type );
        }
        if ( $status != 0 ) {
            die "ERROR: $msg\n";
        }

    }
    else
    {

        ( $pipelineTestXmlPath, $GA_TEST_DATA_PATH ) =
            getPipelineTestXmlpath($GA_TEST_DATA_PATH);
        
        die "Pipeline test XML file is $pipelineTestXmlPath";
        #print Dumper($testsConfRef);
        
        my $testsConfRef = readPipelineTestXml($pipelineTestXmlPath);
        
        my ( $testRef, $testSetRef ) =
            processPipelineTestXml( %$testsConfRef, $testId, $testSetId, $type );
        
        if ( $testRef->{enabled} eq "N" || $testSetRef->{enabled} eq "N" ) {
            print "test $testId for $type $testSetId disabled\n";
            return 0;
        }
        
        my $script = $testRef->{script};
        
        if ( $testRef->{setUp} eq "Y" ) {
            setUpTest( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                       $testId, $testSetId, $type );
        }
        
        if ( $script eq "validate" ) {
            my $configuration =
                generateValidateConfig( $GA_TEST_DATA_PATH, $testId, $testSetId,
                                        $type, %$testRef, %$testSetRef );
            $status =
                runValidate( %$testsConfRef, $GA_TEST_DATA_PATH, $testId, $testSetId,
                             $configuration );
        }
        elsif ( $script eq "goat" ) {
            ( $status, $msg ) =
                runGoat( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                         $testId, $testSetId, $type, %$testRef, %$testSetRef );
        }
        elsif ( $script eq "bustard" ) {
            ( $status, $msg ) =
                runBustard( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                            $testId, $testSetId, $type, %$testRef, %$testSetRef );
        }
        elsif ( $script eq "GERALD" ) {
            ( $status, $msg ) =
                runGerald( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                           $testId, $testSetId, $type, %$testRef, %$testSetRef );
        }
        else {
            die "ERROR: unsupported script in [$testId] test\n";
        }
        
        if ( $testRef->{tearDown} eq "Y" ) {
            tearDownTest( %$testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH,
                          $testId, $testSetId, $type );
        }
        if ( $status != 0 ) {
            die "ERROR: $msg\n";
        }
    }

    return $status;
}    # runTest

1;   # says use was ok
__END__

