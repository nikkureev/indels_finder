import matplotlib.pyplot as plt

path1 = r'C:\Theileria\Parva_Annulata.csv'
path2 = r'C:\Theileria\Parva_Orientalis.csv'
path3 = r'C:\Theileria\Annulata_Orientalis.csv'
path4 = r'C:\Theileria\Parva_Annulata_selected.txt'
path5 = r'C:\Theileria\Parva_Orientalis_selected.txt'
path6 = r'C:\Theileria\Annulata_Orientalis_selected.txt'


def grad_count(input, output, level):
    u, v, sum, u_1, v_1, u_len, v_len = 0, 0, 0, 0, 0, 0, 0
    with open(input, 'r') as table:
        with open(output, 'w') as f:
            for lines in table:
                if float(lines.split(',')[2]) > int(level):
                    f.write(lines)
                    u += float(lines.split(',')[11])
                    u_len += int(lines.split(',')[3])
                    sum += 1
                    u_1 += 1
                else:
                    v += float(lines.split(',')[11])
                    v_len += int(lines.split(',')[3])
                    sum += 1
                    v_1 += 1
    print('+++++++++++++++++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++++++++++++++++')
    print('selected: ', end='')
    print(u/u_1, end=' ; ')
    print(u_1, end=' ; ')
    print(u_len/u_1)
    print('non selected: ', end='')
    print(v/v_1, end=' ; ')
    print(v_1, end=' ; ')
    print(v_len /v_1)
#print('% of all length: ', end='')
#print(u/v * 100)


def plot_func(input, e=1):
    y_list = []
    if e == 1:
        for i in range(50, 100):
            sum_in, score_in = 1, 0
            with open(input, 'r') as table:
                for lines in table:
                    if float(lines.split(',')[2]) > int(i):
                        sum_in += 1
                        score_in += float(lines.split(',')[3])
            y_list.append(score_in/sum_in)
    else:
        for i in range(100, 50, -1):
            sum_in, score_in = 1, 0
            with open(input, 'r') as table:
                for lines in table:
                    if float(lines.split(',')[2]) < int(i):
                        sum_in += 1
                        score_in += float(lines.split(',')[3])
            y_list.append(score_in/sum_in)

    return y_list



def KeyFunc(item):
    return item[1]


def graf(file, a):
    list_1 = []
    list_2 = []
    list_3 = []
    list_4 = []
    len_list = []
    with open(file, 'r') as f:
        for lines in f:
            if 65 <= float(lines.split(',')[a]) < 70:
                list_1.append(lines.split(',')[a])
            if 70 <= float(lines.split(',')[a]) < 75:
                list_1.append(lines.split(',')[a])
            if 75 <= float(lines.split(',')[a]) < 80:
                list_1.append(lines.split(',')[a])
            if 80 <= float(lines.split(',')[a]) < 85:
                list_1.append(lines.split(',')[a])
            if 85 <= float(lines.split(',')[a]) < 90:
                list_2.append(lines.split(',')[a])
            if 90 <= float(lines.split(',')[a]) < 95:
                list_3.append(lines.split(',')[a])
            if float(lines.split(',')[a]) >= 95:
                list_4.append(lines.split(',')[a])
    len_list.append(len(list_1))
    len_list.append(len(list_2))
    len_list.append(len(list_3))
    len_list.append(len(list_4))
    fig = plt.figure()
    plt.bar(range(len(len_list)), len_list)
    plt.show()


def gomology_list_maker(file, a, b):
    with open(file, 'r') as f:
        genome_list = []
        number = 0
        for lines in f:
            if not lines.startswith('#'):
                gom_frg = []
                gom_frg.append(number)
                coord = []
                coord.append(int(lines.split(',')[a]))
                coord.append(int(lines.split(',')[b]))
                gom_frg.append(coord)
                genome_list.append(gom_frg)
                number += 1
    #print(genome_list)
    return genome_list


def space_matrix_maker(list, file):
    space_matrix = []
    for i in range(len(list)):
        spacing = []
        for j in range(len(list)):
            l1 = abs(list[i][1][1] - list[j][1][0])
            l2 = abs(list[i][1][1] - list[j][1][1])
            l3 = abs(list[i][1][0] - list[j][1][0])
            l4 = abs(list[i][1][0] - list[j][1][1])
            spacing.append(min(l1, l2, l3, l4))
        space_matrix.append(spacing)
    with open(file, 'w') as f:
        for i in space_matrix:
            for j in i:
                f.write(str(j) + ' ')
            f.write('\n')
    return space_matrix


def visual(matrix, list, out_list):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] < 2000 and matrix[i][j] > 200:
                if abs(list[i][1][1] - list[i][1][0]) > 100 and \
                        abs(list[j][1][1] - list[j][1][0]) > 100:
                    #print(matrix[i][j])
                    #print(i, end=' ')
                    #print(list[i], end=' ')
                    #print(abs(list[i][1][1] - list[i][1][0]), end='\n')
                    #print(j, end=' ')
                    #print(list[j], end=' ')
                    #print(abs(list[j][1][1] - list[j][1][0]), end='\n')
                    #print('--------------')
                    out_list.append(list[i][1])


def list_analysis(list1, list2):
    for i in list1:
        for j in list2:
            if not((i[0] < j[0] and i[1] < j[0] and i[0] < j[1] and i[1] < j[1])\
            or (i[0] > j[0] and i[1] > j[0] and i[0] > j[1] and i[1] > j[1])):
                print(i, end='^^^^^^^')
                print(j)




Annulata_1 = []
Parva_1 = []
Orientalis_1 = []
Annulata_2 = []
Parva_2 = []
Orientalis_2 = []
x_1_list = [i for i in range(50, 100)]
x_2_list = [i for i in range(100, 50, -1)]

grad_count(path1, path4, 80)

plot_func(path1)
f11 = plt.figure()
plt.plot(x_1_list, plot_func(path1))

plot_func(path1)
f12 = plt.figure()
plt.plot(x_2_list, plot_func(path1, 0))

list_3 = gomology_list_maker(path4, 6, 7)
list_6 = gomology_list_maker(path4, 8, 9)
MATRIX3 = space_matrix_maker(list_3, 'C:/Theileria/DEL_FILE_1.txt')
MATRIX4 = space_matrix_maker(list_6, 'C:/Theileria/DEL_FILE_1.txt')
visual(MATRIX3, list_3, Annulata_1)
print('=============================================')
visual(MATRIX4, list_6, Orientalis_1)

grad_count(path2, path5, 70)

plot_func(path2)
f21 = plt.figure()
plt.plot(x_1_list, plot_func(path2))

plot_func(path2)
f22 = plt.figure()
plt.plot(x_2_list, plot_func(path2, 0))

list_2 = gomology_list_maker(path5, 6, 7)
list_5 = gomology_list_maker(path5, 8, 9)
MATRIX1 = space_matrix_maker(list_2, 'C:/Theileria/DEL_FILE_2.txt')
MATRIX2 = space_matrix_maker(list_5, 'C:/Theileria/DEL_FILE_2.txt')
visual(MATRIX1, list_2, Parva_1)
print('=============================================')
visual(MATRIX2, list_5, Orientalis_2)

grad_count(path3, path6, 70)

plot_func(path3)
f31 = plt.figure()
plt.plot(x_1_list, plot_func(path3))

plot_func(path3)
f32 = plt.figure()
plt.plot(x_2_list, plot_func(path3, 0))

list_1 = gomology_list_maker(path6, 6, 7)
list_4 = gomology_list_maker(path6, 8, 9)
MATRIX5 = space_matrix_maker(list_1, 'C:/Theileria/DEL_FILE_3.txt')
MATRIX6 = space_matrix_maker(list_4, 'C:/Theileria/DEL_FILE_3.txt')
visual(MATRIX5, list_1, Parva_2)
print('=============================================')
visual(MATRIX6, list_4, Annulata_2)
#list_3.sort(key=KeyFunc)
#list_6.sort(key=KeyFunc)

list_analysis(Annulata_1, Annulata_2)



#graf('C:\Theileria\Parva_Annulata.csv', 2)

plt.show()
