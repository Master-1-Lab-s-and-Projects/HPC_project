import sys

def my_func(f_1):
    lines_1 = []
    with open(f_1, 'r') as f:
        for l in f.readlines():
            l = l.strip()
            lines_1.append([int(x) for x in l[1:-1].split(',')])
    return lines_1

if __name__ == '__main__':
    f_1 = sys.argv[1]
    f_2 = sys.argv[2]

    lines_1 = my_func(f_1)
    lines_2 = my_func(f_2)
    for l1,l2 in zip(lines_1, lines_2):
        if l1 != l2:
            print("Erreur: ", l1, l2)


