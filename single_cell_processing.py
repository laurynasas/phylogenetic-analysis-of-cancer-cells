def count_differences(str1, str2):
    count = 0
    str1 = str1.split(",")
    str2 = str2.split(",")
    # print str1
    # print str2
    for el1, el2 in zip(str1,str2):
        if el1 !=el2:
            count +=1
    return count

def calc_distance_matrix_from_file(dir):
    sample_data_file = open(dir, 'rw+')
    all_lines = []
    for line in sample_data_file:
        # if len(line)>2:
        all_lines +=[line]
    sample_data_file.close()

    matrix = []
    for i in all_lines:
        distance_row = []
        for j in all_lines:
            distance_row += [count_differences(i,j)]
        matrix += [distance_row]
    return matrix

def plot_2D_similarity_matrix(matrix):
    from matplotlib import mpl, pyplot
    import numpy as np

    # make values from -5 to 5, for this example
    # zvals = np.random.rand(100, 100) * 10 - 5

    # make a color map of fixed colors
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                         ['blue', 'black', 'red'],
                                                         256)
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    img2 = pyplot.imshow(matrix, interpolation='nearest',
                         cmap=cmap2,
                         origin='lower')

    # make a color bar
    pyplot.colorbar(img2, cmap=cmap2)
    fig = pyplot.figure(1)

    fig.savefig("image2.png")
    pyplot.show()

if __name__ == "__main__":
    dir = "/home/laurynas/workspace/individual_project/data/hou/clustered_data/Bernouli_mixture_model_2017-02-06 16:06:42.239805.txt"
    matrix = calc_distance_matrix_from_file(dir)
    plot_2D_similarity_matrix(matrix)
    # target = open("./data/hou/clustered_data/Bernouli_mixture_model_2017-02-06 16:06:42.239805_matrix.txt", 'w+')
    # for el in matrix:
    #     line = ' | '.join(str(x).ljust(2) for x in el) + "\n"
    #     target.write(line)
    # target.close()

