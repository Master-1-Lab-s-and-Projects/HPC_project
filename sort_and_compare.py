import sys

if __name__=="__main__":
    my_file = sys.argv[1]
    initial_arr = [*range(0,12)]
    with open(my_file, 'r') as f:
        for i,s in enumerate(f.readlines()):
            s = s.strip()
            arr = sorted([int(x) for x in s[1:-1].split(',')])
            if arr != initial_arr:
                print("Different at line", str(i))

