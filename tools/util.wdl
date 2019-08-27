## Description: 
##
## This WDL tool includes custom scripts to parse and reformat 
## task output as part of a workflow. 


task get_read_groups {
   File bam
 
   command {
        samtools view -H ${bam} | grep "@RG" > stdout.txt
   }

   runtime {
       docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
   }

   output { 
       String out = read_string("stdout.txt")
   }
}

task prepare_read_groups_for_star {
   String read_groups

   command <<< 
       echo "${read_groups}" | cut -f 2- | sed -e 's/\t/ /g' | awk '{print}' ORS=' , '| sed 's/ , $//' > stdout.txt
   >>>

   runtime {
       docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
   }

   output { 
       String out = read_string("stdout.txt")
   }
}
