__author__ = 'srio'

import json


class GFile(object):

    def __init__(self):

        self.source = None
        self.oe_list = []


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

    def load_source(self,filename="start.00"):

        self.source = self.load_gfile(filename)

    def load_oe(self,filename="start.01",append=True):
        if append:
            self.oe_list.append( self.load_gfile(filename) )
        else:
            self.oe_list = [self.load_gfile(filename)]


    @classmethod
    def _dictionary2object(cls,dict1):

        class K(object):
            pass

        k = K()
        for key in dict1.keys():
            setattr(k, key, dict1[key])
        return k

    def get_source_as_dictionary(self):
        return self.source

    def get_oe_as_dictionary(self,oe_index=0):
        return self.oe_list[oe_index]

    def get_source_as_object(self):
        return self._dictionary2object(self.get_source_as_dictionary())

    def get_oe_as_object(self,oe_index=0):
        return self._dictionary2object(self.get_oe_as_dictionary(oe_index))

    def get_source_as_json(self):
        return json.dumps(self.get_source_as_dictionary(), sort_keys=True, indent=4)

    def get_oe_as_json(self,oe_index=0):
        return json.dumps(self.get_oe_as_dictionary(oe_index), sort_keys=True, indent=4)

if __name__ == "__main__":

    g = GFile()

    g.load_source(filename="start.00")

    print("\nDICT: \n",g.get_source_as_dictionary())

    print("\nOBJECT: \n",dir(g.get_source_as_object()))

    print("\nJSON: \n",g.get_source_as_json())


    # g.load_oe(filename="start.01")

    # print(g.oe_list[0])
    #
    # # for k in g.oe[0].keys():
    # #     print(k," = ",g.oe[0][k])
    # #
    # # print(json.dumps(g.source, sort_keys=True, indent=4))
    #
    # print(json.dumps(g.oe_list[0], sort_keys=True, indent=4))
    #
    #
    # a = g.get_source_as_object()
    # # print(dir(a))

