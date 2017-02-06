


if __name__ == '__main__':
    aln = AlignIO.read('./data/genotype.phy', 'phylip')
    print aln
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    print dm
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    print(tree)
    draw_ascii(tree)