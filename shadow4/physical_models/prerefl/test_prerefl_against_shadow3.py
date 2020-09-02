import os
from shadow4.physical_models.prerefl.prerefl import PreRefl



if __name__ == "__main__":
    #
    # shadow3 fortran
    #

    # from shadow4.compatibility.global_definitions import SHADOW3_BINARY
    # f = open("prerefl_test.inp",'w')
    # f.write("prerefl_test\nBe5_30.dat\n10000")
    # f.close()
    # print("File written to disk: prerefl_test.inp")
    # os.system("%s < prerefl_test.inp" % SHADOW3_BINARY)


    #
    # python implementation
    #


    prerefl_file = "Be5_30.dat"
    PreRefl.prerefl(interactive=False, SYMBOL="Be",DENSITY=1.848,FILE=prerefl_file,E_MIN=5000.0,E_MAX=30000.0,E_STEP=100.0)


    a = PreRefl()
    a.read_preprocessor_file(prerefl_file)

    refraction_index = a.get_refraction_index(10000.0,verbose=True)