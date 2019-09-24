import os
from shadow4.physical_models.prerefl.prerefl import PreRefl

prerefl_test_dot_inp = """prerefl_test
Be5_30.dat
10000
"""

if __name__ == "__main__":
    #
    # shadow3
    #

    f = open("prerefl_test.inp",'w')
    f.write(prerefl_test_dot_inp)
    f.close()
    print("File written to disk: prerefl_test.inp")

    os.system("/home/manuel/OASYS1.2/shadow3/shadow3 < prerefl_test.inp")


    #
    # python implementation
    #

    a = PreRefl()

    prerefl_file = "Be5_30.dat"
    a.read_preprocessor_file(prerefl_file)


    refraction_index = a.get_refraction_index(10000.0,verbose=True)