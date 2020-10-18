ORFS_FOLDER = 'orfs/'


def read_orffinder(name):
    res = []
    with open(name) as f:
        for line in f:
            s_line = line.split()
            if s_line[0].isnumeric():
                a, b = map(int, s_line[:2])
                if a > b:
                    res.append((b, a, '-'))
                else:
                    res.append((a, b, '+'))
    return set(res)


def read_gff(name):
    res = []
    with open(name) as f:
        for line in f:
            s_line = line.split('\t')
            if len(s_line) > 6 and s_line[3].isnumeric() and s_line[4].isnumeric():
                data = line.split('\t')
                a, b, dir = int(data[3]), int(data[4]), data[6]
                res.append((a, b, dir))
    return set(res)


def imitate_gff(name, segments):
    with open(name + '.gff', 'w') as f:
        f.write('##gff-version\t3\n')
        for num, orf in enumerate(sorted(segments)):
            f.write('Tepidisphaera\tProdigal_v2.6.3\tCDS\t' + str(orf[0]) + '\t' +
                    str(orf[1]) + '\t0\t' + orf[2] + '\t0\tID=1_' + str(num + 1) +
                    ';start=' + str(orf[0]) + ';end=' + str(orf[1]) + ';\n')


def get_segment(name, segment):
    with open(name) as f:
        sequence = ''.join(f.read().split('\n')[1:])
        return sequence[segment[0] - 1: segment[1]]


def delete_intersect(segments, inter):
    res = set(segments)
    for orf in segments:
        for intt in inter:
            if intt[0] <= orf[0] <= intt[1] or intt[0] <= orf[1] <= intt[1]:
                res.remove(orf)
                break
    return res


def find_operons(segments: set):
    segments = sorted(segments)
    result = []
    for i in range(len(segments) - 1):
        j = i + 1
        if segments[i][2] == segments[j][2] and segments[i][1] + 150 >= segments[j][0]:
            segments[j] = (segments[i][0], segments[j][1], segments[j][2])
        else:
            result.append(segments[i])
    result.append(segments[-1])
    return set(result)


if __name__ == '__main__':
    orffinder = read_orffinder(ORFS_FOLDER + 'ORFinder.ft')
    print('ORFinder', len(orffinder), sorted(orffinder))
    prodigal = read_gff(ORFS_FOLDER + 'Prodigal.gff')
    print('Prodigal', len(prodigal), sorted(prodigal))
    glimmer = read_gff(ORFS_FOLDER + 'Glimmer.gff')
    print('Glimmer', len(glimmer), sorted(glimmer))
    genemark = read_gff(ORFS_FOLDER + 'Genemark.gff')
    print('Genemark', len(genemark), sorted(genemark))
    prokka = read_gff(ORFS_FOLDER + 'PROKKA.gff')
    print('PROKKA', len(prokka), sorted(prokka))
    all_inter = prokka.intersection(orffinder.intersection(prodigal.intersection(glimmer.intersection(genemark))))
    print('Perfect intersections', len(all_inter), all_inter)

    imitate_gff('ORFinder', orffinder)
    imitate_gff('Glimmer', glimmer)

    imitate_gff('Perfect', all_inter)
    operons = find_operons(prokka)
    imitate_gff('Operons', find_operons(prokka))
    print('Operons', len(operons), sorted(operons))
