# simpy.py
import python.pyapi as pyapi
from python import plt_kwargs


if __name__ == "__main__":
    import sys
    import time

    # Input and output file names
    try:
        input_file = sys.argv[1]
        try:
            output_file = sys.argv[2]
        except IndexError:
            output_file = "./output.txt"
    except IndexError:
        input_file = "./input.txt"
        output_file = "./output.txt"

    start = time.time()
    spec = pyapi.calc(input_file=input_file, output_file=output_file)
    print 'Calculation finished.'
    print "Time used:", time.time() - start, 'sec.'
    spec.show(**plt_kwargs)


    # import matplotlib.pyplot as plt
    # plt.plot(spec.y_sep.T)
    # plt.plot(spec.y_vec)
    # plt.yscale('log')
    # plt.ylim(1e-14,1e-4)
    # plt.legend(spec.labels+['total'])
    # plt.show()