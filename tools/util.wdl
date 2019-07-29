task cat { 
    File infile
    
    command {
        cat ${infile}
    }
   
    output { 
        String out = read_string(stdout())
    }
}

task get_read_groups {
   File bam
 
   command {
        samtools view -H ${bam} | grep "@RG"
   }
   output { 
       String out = read_string(stdout())
   }
}

task prepare_read_groups_for_star {
   String read_groups

   command <<< 
       echo "${read_groups}" | cut -f 2- | sed -e 's/\t/ /g' | awk '{print}' ORS=' , '| sed 's/ , $//'
   >>>

   output { 
       String out = read_string(stdout())
   }
}
