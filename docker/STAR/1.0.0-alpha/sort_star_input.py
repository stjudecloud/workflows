import argparse
import os


def sort_lists(target_list, key_list):
    for i in range(1, len(key_list)):
        key = key_list[i]
        target = target_list[i]
        j = i - 1
        while j >= 0 and key < key_list[j]:
            key_list[j+1] = key_list[j]
            target_list[j+1] = target_list[j]
            j -= 1
        key_list[j+1] = key
        target_list[j+1] = target

def sort_fastqs(read_one_fastqs, read_two_fastqs):
    read_one_basenames = [fastq.split(os.sep)[-1] for fastq in read_one_fastqs]
    sort_lists(read_one_fastqs, read_one_basenames)

    read_two_basenames = [fastq.split(os.sep)[-1] for fastq in read_two_fastqs]
    sort_lists(read_two_fastqs, read_two_basenames)

    return (read_one_fastqs, read_two_fastqs)

def sort_read_groups(read_groups_string):
    read_groups = [rg for rg in read_groups_string.split(' , ')]
    rgids = [flags.split(' ')[0].split(':')[1] for flags in read_groups]

    sort_lists(read_groups, rgids)

    return (read_groups, rgids)

def validate(read_one_fastqs, read_two_fastqs, rgids):
    if (len(read_one_fastqs) != len(rgids)):
        raise argparse.ArgumentError(
            'Must have same number of read groups as fastq pairs'
        )
    
    for i, id in enumerate(rgids):
        if (id not in read_one_fastqs[i]) or (id not in read_two_fastqs[i]):
            raise SystemExit(
                'Error: There\'s a mismatch between '
                'read group IDs and fastq file names'
            )

def write_outfiles(read_one_fastqs, read_two_fastqs, read_groups):
    read_one_file = open('read_one_fastqs_sorted.txt', 'w')
    read_one_file.write(','.join(read_one_fastqs))
    read_one_file.close()

    read_two_file = open('read_two_fastqs_sorted.txt', 'w')
    read_two_file.write(','.join(read_two_fastqs))
    read_two_file.close()

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
    validate(read_one_fastqs, read_two_fastqs, rgids)
    write_outfiles(read_one_fastqs, read_two_fastqs, read_groups)