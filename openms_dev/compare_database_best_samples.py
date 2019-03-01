import csv
from itertools import product

def float_eq(f1, f2, delta):
    if abs(f1 - f2) < delta:
        return True

    return False

if __name__ == "__main__":
    fn_main = "mz_database.txt"
    fn_search = "A10_pos_best_samples.txt"

    with open(fn_main,'r') as main_f:
        data_main = csv.reader(main_f, delimiter = " ")
        with open(fn_search,'r') as search_f:
            data_search = csv.reader(search_f, delimiter = " ")
            delta = 0.01

            for m, s in product(data_main, data_search):
                    if float_eq(float(m[2]), float(s[1]), delta):
                        print(m[0], m[1], m[2], s[1], s[0])

