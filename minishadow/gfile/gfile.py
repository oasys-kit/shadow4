__author__ = 'srio'

import json


class GFile(object):

    def __init__(self):

        self.start00 = None
        self.oe = []


    @classmethod
    def load_gfile(self,filename="start.00"):


        fp = open(filename) # Open file on read mode
        lines = fp.read().split("\n") # Create a list containing all lines
        fp.close() # Close file

        out_dictionary = {}
        for line in lines:
            mylist = line.split("=")
            # print(">>>",len(mylist),mylist)
            if len(mylist) == 2:
                try:
                    value = mylist[1].strip()
                    if "." in value:
                        out_dictionary[mylist[0].strip()] = float(value)
                    else:
                        out_dictionary[mylist[0].strip()] = int(value)
                except:
                    out_dictionary[mylist[0].strip()] = (mylist[1].strip())


        return out_dictionary

    def load_start00(self,filename="start.00"):

        self.start00 = self.load_gfile(filename)

    def load_start01(self,filename="start.01"):
        self.oe.append( self.load_gfile(filename) )


if __name__ == "__main__":

    g = GFile()

    g.load_start00(filename="start.00")

    g.load_start01(filename="start.01")

    print(g.oe[0])

    # for k in g.oe[0].keys():
    #     print(k," = ",g.oe[0][k])
    #
    # print(json.dumps(g.start00, sort_keys=True, indent=4))

    print(json.dumps(g.oe[0], sort_keys=True, indent=4))


