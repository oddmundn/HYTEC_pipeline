#!/bin/sh
#############################################################################
# HYTEC-pipeline: A shell script that performs pre-processing, tag extraction
# and SSCS generation in the HYTEC method
# Based on samtools, bwa, TagXtractor and SSCScreator
#############################################################################


help()
{
    echo "HYTEC-pipeline.sh: \nA script that calls tag_to_header_ion and SSCScreator_ion in a pipeline"
    echo "\nUsage: HYTEC-pipeline.sh <options> \nOptions:\n -i <infile> \n -o <outfile> \n -t <tag_to_header options (in quotes) \n -c <SSCScreator options (in quotes)"
}

if [ $# -eq 0 ]
then
    help
    exit
fi

#PROCESS COMMAND LINE OPTIONS

while [ $# -gt 0 ]
do
key="$1"

  case $key in
    -h|--help)
	help
	;;    
    -i|--infile)
	INFILE=$2
	shift
	;;
    -o|--outfile)
	OUTFILE=$2
	shift
	;;
    -t|--target-to-header-options)
	TagXtractor_OPTIONS=$2
	shift # past argument
	;;
    -c|--SSCScreator-options)
    SSCScreator_OPTIONS="$2"
    shift # past argument
    ;;
    *)
    # unknown option
	echo "Unknown options: $1\n"
	help
    ;;
  esac
shift # past argument or value
done

#echo "Infile: $INFILE"
#echo "TagXtractor options: $TagXtractor_OPTIONS"

if ! [ -f $INFILE ]
then
    echo "Input file does not exist"
    exit
fi

FWDOUTPREFIX="${INFILE}.fwd.tagx"
REVOUTPREFIX="${INFILE}.rev.tagx"
FWDREPORTFILE="${INFILE}.fwd.SSCS.report"
REVREPORTFILE="${INFILE}.rev.SSCS.report"

#Pipeline is run for forward strand reads frist
echo "#HYTEC-PIPELINE: Forward strand TAG_TO_HEADER, Alignment and CONSENSUSMAKER is running#"
samtools view -b -F20 $INFILE | #-F20 means all mapped reads not reverse (4+16)
samtools fastq - |
    TagXtractor.py --infile - --outprefix $FWDOUTPREFIX --taglen 6 --spacerlen 14  --filtspacer AGCTCTTCCGATCT $TagXtractor_OPTIONS |
    bwa mem -t 4 /mnt/SEQ/ref/hg19.fasta - |
    samtools view -bu - | samtools sort -T "${INFILE}" - |
    SSCScreator.py --infile - --outfile "${INFILE}.fwd.fq.tmp" --tagprefix "${INFILE}.fwd" --reportfile $FWDREPORTFILE  --filt sn --remove_false "${INFILE}.fwd.falsefam" --Ncutoff 0.3 $SSCScreator_OPTIONS

#Then run for reverse strand
echo "#HYTEC-PIPELINE: Reverse strand TagXtractor, Alignment and CONSENSUSMAKER is running#"
samtools view -b -f16 $INFILE |    # -f means only reverse reads
samtools fastq - |
    TagXtractor.py --infile - --outprefix $REVOUTPREFIX --taglen 6 --spacerlen 14 --filtspacer AGCTCTTCCGATCT $TagXtractor_OPTIONS |
    bwa mem -t 4 /mnt/SEQ/ref/hg19-rev.fasta - |     #Use reverse reference gen.
    samtools view -bu - | samtools sort -T "${INFILE}" - |
    SSCScreator.py --infile - --outfile "${INFILE}.rev.fq.tmp" --tagprefix "${INFILE}.rev" --reportfile $REVREPORTFILE  --filt sn --remove_false "${INFILE}.rev.falsefam" --Ncutoff 0.3 $SSCScreator_OPTIONS

#Concatenate the SSCS fq files
echo "#HYTEC-PIPELINE: Concatenating SSCS fastq files#"

#Check that SSCS tmp files exist
if ! [ -f ${INFILE}.fwd.fq.tmp ]
then
    echo "Fwd SSCS fastq file did not exist"
    rm ${INFILE}.rev.fq.tmp
    exit
fi

if ! [ -f ${INFILE}.rev.fq.tmp ]
then
    echo "Reverse SSCS fastq file did not exist"
        rm ${INFILE}.fwd.fq.tmp
    exit
fi

cat ${INFILE}.fwd.fq.tmp ${INFILE}.rev.fq.tmp > ${INFILE}.fq.tmp
rm ${INFILE}.fwd.fq.tmp
rm ${INFILE}.rev.fq.tmp

#Align and sort
echo "#HYTEC-PIPELINE: Final aligment and sorting#"
bwa mem /mnt/SEQ/ref/hg19.fasta ${INFILE}.fq.tmp |
    samtools view -bu - | samtools sort -T "${INFILE}" - > $OUTFILE
rm ${INFILE}.fq.tmp
