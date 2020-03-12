import argparse


def sort_fastqs(read_one_fastqs, read_two_fastqs):
    read1_names = [fastq.split('/')[-1] for fastq in read_one_fastqs]
    read2_names = [fastq.split('/')[-1] for fastq in read_two_fastqs]
    
    for i in range(1, len(read1_names)):
        key = read1_names[i]
        key_pair = read_one_fastqs[i]
        j = i - 1
        while j >= 0 and key < read1_names[j]:
            read1_names[j+1] = read1_names[j]
            read_one_fastqs[j+1] = read_one_fastqs[j]
            j -= 1
        read1_names[j+1] = key
        read_one_fastqs[j+1] = key_pair
    
    for i in range(1, len(read2_names)):
        key = read2_names[i]
        key_pair = read_two_fastqs[i]
        j = i - 1
        while j >= 0 and key < read2_names[j]:
            read2_names[j+1] = read2_names[j]
            read_two_fastqs[j+1] = read_two_fastqs[j]
            j -= 1
        read2_names[j+1] = key
        read_two_fastqs[j+1] = key_pair
      
    return (read_one_fastqs, read_two_fastqs)

def sort_read_groups(read_groups_string):
    read_groups = [rg for rg in read_groups_string.split(' , ')]
    rgids = [flags.split(' ')[0].split(':')[1] for flags in read_groups]

    for i in range(1, len(rgids)):
        key = rgids[i]
        key_pair = read_groups[i]
        j = i - 1
        while j >= 0 and key < rgids[j]:
            rgids[j+1] = rgids[j]
            read_groups[j+1] = read_groups[j]
            j -= 1
        rgids[j+1] = key
        read_groups[j+1] = key_pair

    return (read_groups, rgids)

def validate(read_one_fastqs, read_two_fastqs, rgids):
    if (len(read_one_fastqs) != len(rgids)):
        raise argparse.ArgumentError(
            'Must have same number of read groups as fastq pairs'
        )
    
    for i, id in enumerate(rgids):
        if (id not in read_one_fastqs[i]) or (id not in read_two_fastqs[i]):
            raise SystemExit('Error: read group id not in fastqs')

def write_outfiles(read_one_fastqs, read_two_fastqs, read_groups):
    read1_file = open('read_one_fastqs_sorted.txt', 'w')
    read1_file.write(','.join(read_one_fastqs))
    read1_file.close()

    read2_file = open('read_two_fastqs_sorted.txt', 'w')
    read2_file.write(','.join(read_two_fastqs))
    read2_file.close()

    read_group_file = open('read_groups_sorted.txt', 'w')
    read_group_file.write(' , '.join(read_groups))
    read_group_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_one_fastqs', nargs='+', required=True)
    parser.add_argument('--read_two_fastqs', nargs='+', required=True)
    parser.add_argument('--read_groups', required=True)

    args = parser.parse_args()
    
    if (len(args.read_one_fastqs) != len(args.read_two_fastqs)):
        raise argparse.ArgumentError(
            'Must have the same number of read one fastqs as read two fastqs'
        )
    read_one_fastqs, read_two_fastqs = sort_fastqs(
        args.read_one_fastqs, args.read_two_fastqs
    )
    read_groups, rgids = sort_read_groups(args.read_groups)
    write_outfiles(read_one_fastqs, read_two_fastqs, read_groups)