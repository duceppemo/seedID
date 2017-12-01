#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
export baseDir=""${HOME}"/analyses/seedID/170822"

#reads
export fastq="/media/6tb_raid10/data/pollen/Pollen_100pM.sample/pollen_nomatch.fastq.gz"

#barcodes
barcodes="/media/6tb_raid10/data/pollen/pollen_barcodes.csv"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=4

#program location
export prog=""${HOME}"/prog"
export picard=""${prog}"/picard-tools/picard.jar"


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"
export memCdHit=$((mem*1000))


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
export logs=""${baseDir}"/logs"
export demul=""${baseDir}"/demul"
export trimmed=""${baseDir}"/trimmed"
export clustered=""${baseDir}"/clustered"
export blast=""${baseDir}"/blast"
export aligned=""${baseDir}"/aligned"
export centrifuge=""${baseDir}"/centrifuge"

#create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$demul" ] || mkdir -p "$demul"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$clustered" ] || mkdir -p "$clustered"
[ -d "$blast" ] || mkdir -p "$blast"
[ -d "$aligned" ] || mkdir -p "$aligned"
[ -d "$centrifuge" ] || mkdir -p "$centrifuge"


######################
#                    #
#   Demultiplexing   #
#                    #
######################


# Process barcode file

function Demultiplex()
{
    # MIRL15-Bnap-01-ITS4_BC81 CCTGCCATTCGCGAT
    name=$(cut -d "," -f 1 <<< "$1")
    seqence=$(cut -d "," -f 2 <<< "$1")

    #find barcodes and trim them
    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in="$fastq" \
        literal="$seqence" \
        k=17 hdist=1 \
        outm="${demul}"/"${name}".fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/"${name}".txt)
}

#make function available to parallel
export -f Demultiplex  # -f is to export functions

cat "$barcodes" \
    | parallel  --env Demultiplex \
                --env cpu \
                --env maxProc \
                --env memJava \
                --env fastq \
                --env demul \
                --env logs \
                --jobs "$maxProc" \
                "Demultiplex {}"

# make fasta with barcode file
cat "$barcodes" | sed -e 's/^/>/' -e 's/,/\n/' > "${baseDir}"/barcodes.fasta

#nomatch
bbduk.sh "$memJava" \
    threads="$cpu" \
    in="$fastq" \
    ref="${baseDir}"/barcodes.fasta \
    k=17 hdist=1 \
    out="${demul}"/nomatch.fastq.gz \
    pigz=t \
    unpigz=t \
    2> >(tee "${logs}"/nomatch.txt)

# Total Removed:              6937501 reads (53.96%)  2153420965 bases (56.68%)
for i in $(find "$logs" -type f -name "*.txt" | grep -v "nomatch"); do
    name=$(basename "${i%.txt}")
    entry_number=$(cat "$i" | grep -F 'Total Removed:' | tr " " "\t" | awk '{print $3}')
    echo -e ""$name"\t"$entry_number"" | tee -a "${logs}"/demul_stats.tsv
done

echo -e "nomatch\t$(cat "${logs}"/nomatch.txt \
    | grep -F "Result" \
    | tr " " "\t" \
    | awk '{print $2}')" \
    | tee -a "${logs}"/demul_stats.tsv

#Sort file
cat "${logs}"/demul_stats.tsv \
    | sort -k1,1 \
    > "${logs}"/tmp

#rename file
mv "${logs}"/tmp "${logs}"/demul_stats.tsv


# Block comment
: <<'END'
#### Sabre


#Make barcode file for sabre
for i in $(cat "$barcodes"); do
    name=$(cut -d "," -f 1 <<< "$i")
    seqence=$(cut -d "," -f 2 <<< "$i")

    echo -e ""$seqence"\t""${demul}"/"${name}".fastq \
        >> "${baseDir}"/barcodes_sabre.txt
done

#sabre
sabre se \
    -f "$fastq" \
    -b "${baseDir}"/barcodes_sabre.txt \
    -u "${demul}"/nomatch.fastq \
    -m 1 | tee "${logs}"/sabre.txt

#replace name
 awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}1' \
    /home/bioinfo/analyses/plant_demul/barcodes_sabre.txt \
    "${logs}"/sabre.txt \
    > "${logs}"/sabre.txt.tmp

#  mv "${logs}"/sabre.txt.tmp "${logs}"/sabre.txt
END


#########################
#                       #
#   Barcode trimming    #
#                       #
#########################


#TODO -> trim minlen here or not?

# Trim barcode sequences
function TrimBarcodes()
{
    sample=$(basename "${1%.fastq.gz}")

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in="$1" \
        ref="/media/6tb_raid10/data/plant/ITS2/reference/barcodes.fasta" \
        ktrim=l k=17 hdist=1 \
        qtrim=lr trimq=10 \
        minlen=64 \
        out="${trimmed}"/"${sample}"_trimmed.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/"${sample}"_trimming.txt)
}

#make function available to parallel
export -f TrimBarcodes  # -f is to export functions

find "$demul" -type f -name "*.fastq.gz" \
    | parallel  --env TrimBarcodes \
                --env cpu \
                --env maxProc \
                --env memJava \
                --env baseDir \
                --env trimmed \
                --env logs \
                --jobs "$maxProc" \
                "TrimBarcodes {}"




#######################
#                     #
#   Read filtering    #
#                     #
#######################


# TODO
# remove reads too long, for sure, maybe the too short ones too???


############################
#                          #
#   Sequence clustering    #
#                          #
############################


#CD-HIT-EST
function Cluster()
{
    sample=$(basename "${1%_trimmed.fastq.gz}")

    #convert fastq to fasta and pipe into cd-hit
    zcat "$1" | sed -n '1~4s/^@/>/p;2~4p' > "${clustered}"/"${sample}".fasta

    #run cd-hitwith 97% identity and 80% length
    cd-hit-est \
        -i "${clustered}"/"${sample}".fasta \
        -o "${clustered}"/"${sample}"_clustered.fasta \
        -c 0.97 \
        -s 0.90 \
        -T $((cpu/maxProc)) \
        -M $((memCdHit/maxProc)) \
        -n 10 \
        -d 0

    # Modify fasta headers to include clustering info
    # filter out the clusters with less than 20 reads
    # min cluster lenght = 200 and max = 400
    perl "${HOME}"/scripts/cd-hitHeaderRenamer.pl \
        "${clustered}"/"${sample}"_clustered.fasta \
        "${clustered}"/"${sample}"_clustered.fasta.clstr \
        "${clustered}"/"${sample}"_clustered_mod.fasta \
        20 200 400

    #sort by most abundant first
    cat "${clustered}"/"${sample}"_clustered_mod.fasta \
        | sed -e '/>/s/^/@/' -e '/>/s/$/#/' \
        | tr -d "\n" \
        | tr "@" "\n" \
        | sort -t "=" -nrk2 \
        | tr "#" "\n" \
        | sed -e '/^$/d' \
        > "${clustered}"/"${sample}"_sorted.fasta
}

# make function available to parallel
export -f Cluster  # -f is to export functions

# Run in parallel on multiple samples
find "$trimmed" -type f -name "*.fastq.gz" \
    | parallel  --env Cluster \
                --env cpu \
                --env maxProc \
                --env memCdHit \
                --env clustered \
                --env logs \
                --jobs "$maxProc" \
                "Cluster {}"


##################
#                #
#     BLAST      #
#                #
##################


# Blast clustered sequeces on nt
function Blast()
{
    sample=$(basename "${1%_sorted.fasta}")

    blastn \
        -db "/media/6tb_raid10/db/plant_ITS/2017-10-16_plant_ITS.fasta" \
        -query "$1" \
        -num_threads $((cpu/maxProc)) \
        -out "${blast}"/"${sample}".blastn \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -max_target_seqs 1

    # Add header to blast output
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        > "${blast}"/"${sample}".blastn.tmp
    cat "${blast}"/"${sample}".blastn >> "${blast}"/"${sample}".blastn.tmp
    mv "${blast}"/"${sample}".blastn.tmp "${blast}"/"${sample}".blastn
}


# make function Cluster to parallel
export -f Blast  # -f is to export functions

# Run in parallel on multiple samples
find "$clustered" -type f -name "*_sorted.fasta" \
    | parallel  --env Blast \
                --env cpu \
                --env maxProc \
                --env blast \
                --jobs "$maxProc" \
                "Blast {}"


# Blast clustered sequeces on regulated species db
function BlastRegulated()
{
    sample=$(basename "${1%_sorted.fasta}")

    blastn \
        -db /media/6tb_raid10/data/plant/ITS2/reference/2017-03-29_invasive-plant-species_ITS.fasta \
        -query "$1" \
        -num_threads $((cpu/maxProc)) \
        -out "${blast}"/"${sample}"_regulated.blastn \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -max_target_seqs 1

    # Add header to blast output
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        > "${blast}"/"${sample}"_regulated.blastn.tmp
    cat "${blast}"/"${sample}"_regulated.blastn >> "${blast}"/"${sample}"_regulated.blastn.tmp
    mv "${blast}"/"${sample}"_regulated.blastn.tmp "${blast}"/"${sample}"_regulated.blastn
}


# make function available to parallel
export -f BlastRegulated  # -f is to export functions

# Run in parallel on multiple samples
find "$clustered" -type f -name "*_sorted.fasta" \
    | parallel  --env BlastRegulated \
                --env cpu \
                --env maxProc \
                --env blast \
                --jobs "$maxProc" \
                "BlastRegulated {}"


###############
#             #
#   Mapping   #
#             #
###############


# Trying to use a different strategy that with the blast

export genome="/media/6tb_raid10/data/plant/ITS2/reference/crops_ITS2_refseq.fasta"

#index reference genome for bwa
if [ -e "${genome}".sa ] && [ -e "${genome}".amb ] \
&& [ -e "${genome}".ann ] && [ -e "${genome}".pac ] \
&& [ -e "${genome}".bwt ]; then #check if indexing already done
    echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing reference genome $(basename "$genome") with bwa..."
    bwa index "$genome"
fi

#index reference genome for samtools
if [ -e "${genome}".fai ]; then #check if indexing already done
    echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing reference genome $(basename "$genome") with Samtools..."
    samtools faidx "$genome"
fi


function Align()
{
    sample=$(basename "${1%_trimmed.fastq.gz}")

    #Align trimmed reads to reference genome (SAM output with Read Group information),
    #convert SAM to BAM, filter out unmapped and sort
    bwa mem -t $((cpu/maxProc)) -r 1 -a -M "$genome" "$1" | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) -m 10G -o "${aligned}"/"${sample}".bam -

    sambamba markdup -r -t $((cpu/maxProc)) "${aligned}"/"${sample}".bam "${aligned}"/"${sample}"_nodup.bam

    # remove duplicates for sinle-end reads
    java -Xmx64g -jar "$picard" MarkDuplicates \
        INPUT="${aligned}"/"${sample}".bam \
        OUTPUT="${aligned}"/"${sample}"_nodup.bam \
        METRICS_FILE="${aligned}"/"${sample}"_duplicates.txt \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true

    # Merge and sort all merged files
    samtools sort -@ $((cpu/maxProc)) -m 10G \
        -o "${aligned}"/"${sample}"_nodup_sorted.bam \
        "${aligned}"/"${sample}"_nodup.bam

    # index bam file
    samtools index "${aligned}"/"${sample}"_nodup_sorted.bam

    find "$aligned" -type f ! -name "*sorted.bam" -exec rm {} \;
}

# make function available to parallel
export -f Align  # -f is to export functions

# Run in parallel on multiple samples
find "$trimmed" -type f -name "*.fastq.gz" \
    | parallel  --env Align \
                --env cpu \
                --env maxProc \
                --env genome \
                --env aligned \
                --jobs "$maxProc" \
                "Align {}"


# TODO -> What is next???
# check if coverage is good over all reference sequence length



##################
#                #
#   Centrifuge   #
#                #
##################


# Trying to use a different strategy that with the blast

# use nt database
# db="/media/6tb_raid10/db/centrifuge/nt"
db="/media/6tb_raid10/db/centrifuge/2017-10-17_plant_ITS"

for i in $(find "$trimmed" -type f -name "*.fastq.gz"); do

    sample=$(basename "${i%_trimmed.fastq.gz}")
    [ -d "${centrifuge}"/"${sample}" ] || mkdir -p "${centrifuge}"/"${sample}"

    centrifuge \
        -q \
        -p $(nproc) \
        -t \
        --seed "$RANDOM" \
        -x "$db" \
        -U "$i" \
        --report-file "${centrifuge}"/"${sample}"/"${sample}"_report.tsv \
        > "${centrifuge}"/"${sample}"/"${sample}".tsv

    #Prepare result for display with Krona
    cat "${centrifuge}"/"${sample}"/"${sample}".tsv | \
        cut -f 1,3 | \
        ktImportTaxonomy /dev/stdin -o  "${centrifuge}"/"${sample}"/"${sample}".html

    #visualize the resutls in Firefow browser
    # firefox file://"${centrifuge}"/"${sample}"/"${sample}".html &

done


# TODO
    # extract reads that match unexpected species and compare them to the reads that match to the expected crop



################
#              #
#   Distance   #
#              #
################



#distance matrix for all invasive species and crops
cat /media/6tb_raid10/data/plant/ITS2/reference/crops_ITS2_refseq.fasta \
    /media/6tb_raid10/data/plant/ITS2/reference/2017-03-29_invasive-plant-species_ITS_cleaned.fasta | \
muscle \
    -in /dev/stdin \
    -out /media/6tb_raid10/data/plant/ITS2/reference/2017-03-29_all_ITS2_alignment.fasta

