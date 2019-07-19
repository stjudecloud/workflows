task infer_experiment {
    File bam
    Int? sample_size
    Int? map_qual
    File refgene_bed 
 
    command {
        infer_experiment.py -i ${bam} -r ${refgene_bed} ${"-s" + sample_size} ${"-q" + map_qual} 
    }

    output {
       String out = read_string(stdout())
    }
}
